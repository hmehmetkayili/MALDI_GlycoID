"""HTTP routing, security middleware, and TestClient compatibility engine for MALDI-GlycoID (Stage 12).

Provides a unified FastAPI/Starlette-compatible API layer (`FastAPI`, `APIRouter`,
`HTTPException`, `Request`, `Response`, `JSONResponse`, `HTMLResponse`, `UploadFile`,
`TestClient`, and `run_http_server`) that operates seamlessly both in environments
with `fastapi`/`uvicorn` installed and in air-gapped scientific Python distributions
using only the Python 3.11 standard library.
"""

from __future__ import annotations

import asyncio
from dataclasses import dataclass, field
from email.parser import BytesParser
from email.policy import default as default_email_policy
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import inspect
import io
import json
_STATIC_MIME_TYPES: Dict[str, str] = {
    ".css": "text/css; charset=utf-8",
    ".js": "application/javascript; charset=utf-8",
    ".html": "text/html; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".svg": "image/svg+xml",
    ".ico": "image/x-icon",
    ".txt": "text/plain; charset=utf-8",
}
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple, Union
from urllib.parse import parse_qs, urlencode, urlsplit

from app.config import settings

SECURITY_HEADERS: Dict[str, str] = {
    "X-Content-Type-Options": "nosniff",
    "X-Frame-Options": "DENY",
    "Referrer-Policy": "strict-origin-when-cross-origin",
    "X-XSS-Protection": "1; mode=block",
    "Content-Security-Policy": (
        "default-src 'self'; "
        "img-src 'self' data:; "
        "style-src 'self' 'unsafe-inline'; "
        "script-src 'self';"
    ),
}


class HTTPException(Exception):
    """HTTP error exception carrying status code and structured detail message."""

    def __init__(self, status_code: int, detail: str) -> None:
        super().__init__(detail)
        self.status_code = int(status_code)
        self.detail = str(detail)


@dataclass
class UploadFile:
    """Uploaded file container supporting both sync and async `.read()`."""

    filename: str
    content: bytes
    content_type: str = "application/octet-stream"

    def read_bytes(self) -> bytes:
        return self.content

    async def read(self) -> bytes:
        return self.content


class FormDataPayload(dict):
    """Multi-value form dictionary supporting `.get()` and `.getlist()`."""

    def __init__(self) -> None:
        super().__init__()
        self._lists: Dict[str, List[Any]] = {}

    def add(self, key: str, value: Any) -> None:
        self._lists.setdefault(key, []).append(value)
        super().__setitem__(key, value)

    def getlist(self, key: str) -> List[Any]:
        return list(self._lists.get(key, []))


@dataclass
class Response:
    """HTTP response container with status code, headers, and byte/text body."""

    content: Union[bytes, str] = b""
    status_code: int = 200
    media_type: str = "text/plain; charset=utf-8"
    headers: Dict[str, str] = field(default_factory=dict)

    @property
    def body(self) -> bytes:
        if isinstance(self.content, bytes):
            return self.content
        return str(self.content).encode("utf-8")

    def finalized_headers(self) -> Dict[str, str]:
        merged = dict(SECURITY_HEADERS)
        for k, v in self.headers.items():
            merged[k] = v
        if not any(k.lower() == "content-type" for k in merged):
            merged["Content-Type"] = self.media_type
        return merged


class JSONResponse(Response):
    """JSON-serialized HTTP response."""

    def __init__(
        self,
        content: Any,
        status_code: int = 200,
        headers: Optional[Dict[str, str]] = None,
    ) -> None:
        serialized = json.dumps(content, ensure_ascii=False).encode("utf-8")
        super().__init__(
            content=serialized,
            status_code=status_code,
            media_type="application/json",
            headers=dict(headers or {}),
        )


class HTMLResponse(Response):
    """HTML HTTP response."""

    def __init__(
        self,
        content: Union[str, bytes],
        status_code: int = 200,
        headers: Optional[Dict[str, str]] = None,
    ) -> None:
        super().__init__(
            content=content,
            status_code=status_code,
            media_type="text/html; charset=utf-8",
            headers=dict(headers or {}),
        )


class Request:
    """Incoming HTTP request representation."""

    def __init__(
        self,
        method: str,
        path: str,
        query_params: Optional[Dict[str, str]] = None,
        headers: Optional[Mapping[str, str]] = None,
        body: bytes = b"",
        preparsed_form: Optional[FormDataPayload] = None,
    ) -> None:
        self.method = method.upper()
        self.path = path
        self.query_params: Dict[str, str] = dict(query_params or {})
        self.headers: Dict[str, str] = {str(k).lower(): str(v) for k, v in (headers or {}).items()}
        self.body = body
        self._preparsed_form = preparsed_form

    def json_sync(self) -> Dict[str, Any]:
        if not self.body:
            return {}
        try:
            data = json.loads(self.body.decode("utf-8"))
            return data if isinstance(data, dict) else {"_value": data}
        except Exception as exc:
            raise HTTPException(status_code=400, detail=f"Invalid JSON request payload: {exc}") from exc

    async def json(self) -> Dict[str, Any]:
        return self.json_sync()

    def form_sync(self) -> FormDataPayload:
        if self._preparsed_form is not None:
            return self._preparsed_form

        form = FormDataPayload()
        ctype = self.headers.get("content-type", "")
        if "multipart/form-data" in ctype:
            raw_msg = (
                f"Content-Type: {ctype}\r\nMIME-Version: 1.0\r\n\r\n".encode("utf-8")
                + self.body
            )
            msg = BytesParser(policy=default_email_policy).parsebytes(raw_msg)
            if msg.is_multipart():
                for part in msg.iter_parts():
                    cd = part.get("Content-Disposition", "")
                    if "form-data" not in cd:
                        continue
                    field_name = part.get_param("name", header="content-disposition")
                    filename = part.get_filename()
                    payload_bytes = part.get_payload(decode=True) or b""
                    if not field_name:
                        continue
                    if filename is not None:
                        form.add(
                            str(field_name),
                            UploadFile(
                                filename=str(filename),
                                content=payload_bytes,
                                content_type=part.get_content_type() or "application/octet-stream",
                            ),
                        )
                    else:
                        form.add(str(field_name), payload_bytes.decode("utf-8", errors="replace"))
        elif "application/x-www-form-urlencoded" in ctype:
            parsed = parse_qs(self.body.decode("utf-8", errors="replace"), keep_blank_values=True)
            for k, vals in parsed.items():
                for v in vals:
                    form.add(k, v)
        self._preparsed_form = form
        return form

    async def form(self) -> FormDataPayload:
        return self.form_sync()


@dataclass
class RouteEntry:
    """Registered HTTP route entry."""

    method: str
    path: str
    endpoint: Callable[..., Any]
    tags: Tuple[str, ...] = field(default_factory=tuple)


class APIRouter:
    """Modular API router supporting prefix and method decorators."""

    def __init__(self, prefix: str = "", tags: Optional[Sequence[str]] = None) -> None:
        self.prefix = prefix.rstrip("/")
        self.tags = tuple(tags or ())
        self.routes: List[RouteEntry] = []

    def add_api_route(self, path: str, endpoint: Callable[..., Any], methods: Sequence[str]) -> None:
        full_path = f"{self.prefix}/{path.lstrip('/')}" if path != "/" else (self.prefix or "/")
        if not full_path.startswith("/"):
            full_path = "/" + full_path
        for m in methods:
            self.routes.append(
                RouteEntry(
                    method=m.upper(),
                    path=full_path,
                    endpoint=endpoint,
                    tags=self.tags,
                )
            )

    def get(self, path: str) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
            self.add_api_route(path, func, methods=["GET"])
            return func

        return decorator

    def post(self, path: str) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
            self.add_api_route(path, func, methods=["POST"])
            return func

        return decorator


class FastAPI:
    """Application router, static file server, and security middleware engine."""

    def __init__(
        self,
        title: str = "MALDI-GlycoID",
        version: str = "1.0.0",
        description: str = "",
    ) -> None:
        self.title = title
        self.version = version
        self.description = description
        self.routes: List[RouteEntry] = []
        self.static_mounts: Dict[str, Path] = {}

    def include_router(self, router: APIRouter, prefix: str = "") -> None:
        extra = prefix.rstrip("/")
        for r in router.routes:
            path = f"{extra}{r.path}" if extra else r.path
            self.routes.append(
                RouteEntry(
                    method=r.method,
                    path=path,
                    endpoint=r.endpoint,
                    tags=r.tags,
                )
            )

    def mount_static(self, url_prefix: str, directory: Union[str, Path]) -> None:
        norm_prefix = "/" + url_prefix.strip("/")
        self.static_mounts[norm_prefix] = Path(directory).resolve()

    def get(self, path: str) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
            self.routes.append(RouteEntry(method="GET", path=path, endpoint=func))
            return func

        return decorator

    def post(self, path: str) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
            self.routes.append(RouteEntry(method="POST", path=path, endpoint=func))
            return func

        return decorator

    def dispatch(self, request: Request) -> Response:
        """Execute a `Request` against registered routes, static mounts, and security checks."""
        # Enforce maximum upload size limit on request body
        if len(request.body) > settings.max_upload_size_bytes:
            return JSONResponse(
                content={
                    "detail": (
                        f"Payload size ({len(request.body)} bytes) exceeds maximum allowed upload limit "
                        f"({settings.max_upload_size_bytes} bytes)."
                    )
                },
                status_code=413,
            )

        effective_method = "GET" if request.method == "HEAD" else request.method

        # Check static file mounts first for GET/HEAD requests
        if effective_method == "GET":
            for prefix, root_dir in self.static_mounts.items():
                if request.path.startswith(prefix + "/"):
                    rel_part = request.path[len(prefix) + 1 :]
                    # Reject path traversal sequences
                    if ".." in rel_part or rel_part.startswith("/") or "\\" in rel_part:
                        return JSONResponse({"detail": "Invalid static file path."}, status_code=400)
                    target = (root_dir / rel_part).resolve()
                    try:
                        target.relative_to(root_dir)
                    except ValueError:
                        return JSONResponse({"detail": "Access denied outside static directory."}, status_code=403)
                    if not target.is_file():
                        return JSONResponse({"detail": f"Static file not found: {rel_part}"}, status_code=404)
                    mime = _STATIC_MIME_TYPES.get(target.suffix.lower(), "application/octet-stream")
                    return Response(
                        content=target.read_bytes(),
                        status_code=200,
                        media_type=mime,
                    )

        # Match registered API or template route
        matched_Routes = [r for r in self.routes if r.path == request.path]
        if not matched_Routes:
            return JSONResponse({"detail": f"Endpoint not found: {request.path}"}, status_code=404)

        route = next((r for r in matched_Routes if r.method == effective_method), None)
        if route is None:
            return JSONResponse(
                {"detail": f"Method {request.method} not allowed for {request.path}"},
                status_code=405,
            )

        try:
            sig = inspect.signature(route.endpoint)
            kwargs: Dict[str, Any] = {}
            if "request" in sig.parameters:
                kwargs["request"] = request

            result = route.endpoint(**kwargs)
            if inspect.isawaitable(result):
                result = asyncio.run(result)

            if isinstance(result, Response):
                return result
            return JSONResponse(content=result, status_code=200)
        except HTTPException as http_exc:
            return JSONResponse(content={"detail": http_exc.detail}, status_code=http_exc.status_code)
        except Exception as exc:
            return JSONResponse(
                content={"detail": f"Internal processing error: {exc}"},
                status_code=500,
            )


class ClientResponse:
    """Response returned by `TestClient` matching `httpx.Response` / `Starlette` interface."""

    def __init__(self, response: Response) -> None:
        self._response = response
        self.status_code: int = response.status_code
        self.headers: Dict[str, str] = {k.lower(): v for k, v in response.finalized_headers().items()}
        self.raw_headers: Dict[str, str] = response.finalized_headers()
        self.content: bytes = response.body

    @property
    def text(self) -> str:
        return self.content.decode("utf-8", errors="replace")

    def json(self) -> Any:
        return json.loads(self.text)


class TestClient:
    """In-memory synchronous test client for `FastAPI` applications."""

    __test__ = False

    def __init__(self, app: FastAPI) -> None:
        self.app = app

    def _split_url(
        self,
        url: str,
        params: Optional[Mapping[str, Any]] = None,
    ) -> Tuple[str, Dict[str, str]]:
        parsed = urlsplit(url)
        q_dict: Dict[str, str] = {}
        if parsed.query:
            for k, vals in parse_qs(parsed.query, keep_blank_values=True).items():
                if vals:
                    q_dict[k] = str(vals[-1])
        if params:
            for k, v in params.items():
                if v is not None:
                    q_dict[str(k)] = str(v)
        return parsed.path or "/", q_dict

    def get(
        self,
        url: str,
        params: Optional[Mapping[str, Any]] = None,
        headers: Optional[Mapping[str, str]] = None,
    ) -> ClientResponse:
        path, q_dict = self._split_url(url, params=params)
        req = Request(method="GET", path=path, query_params=q_dict, headers=headers)
        return ClientResponse(self.app.dispatch(req))

    def post(
        self,
        url: str,
        json: Optional[Any] = None,
        data: Optional[Mapping[str, Any]] = None,
        files: Optional[ Union[Mapping[str, Any], Sequence[Tuple[str, Any]]] ] = None,
        params: Optional[Mapping[str, Any]] = None,
        headers: Optional[Mapping[str, str]] = None,
        content: Optional[bytes] = None,
    ) -> ClientResponse:
        path, q_dict = self._split_url(url, params=params)
        hdrs: Dict[str, str] = {str(k).lower(): str(v) for k, v in (headers or {}).items()}
        body_bytes = content or b""
        preparsed: Optional[FormDataPayload] = None

        if json is not None:
            import json as _json

            body_bytes = _json.dumps(json).encode("utf-8")
            hdrs.setdefault("content-type", "application/json")
        elif files is not None or data is not None:
            preparsed = FormDataPayload()
            if data:
                for k, v in data.items():
                    if isinstance(v, (list, tuple)):
                        for item in v:
                            preparsed.add(str(k), str(item))
                    elif v is not None:
                        preparsed.add(str(k), str(v))
            if files:
                file_items: List[Tuple[str, Any]] = (
                    list(files.items()) if isinstance(files, Mapping) else list(files)
                )
                total_len = 0
                for field_name, file_spec in file_items:
                    if isinstance(file_spec, tuple):
                        fname = str(file_spec[0])
                        raw_f = file_spec[1]
                        ctype = str(file_spec[2]) if len(file_spec) > 2 else "application/octet-stream"
                    else:
                        fname = getattr(file_spec, "name", "upload.dat")
                        raw_f = file_spec
                        ctype = "application/octet-stream"

                    if isinstance(raw_f, (bytes, bytearray)):
                        f_bytes = bytes(raw_f)
                    elif isinstance(raw_f, str):
                        f_bytes = raw_f.encode("utf-8")
                    elif hasattr(raw_f, "read"):
                        read_val = raw_f.read()
                        f_bytes = read_val.encode("utf-8") if isinstance(read_val, str) else bytes(read_val)
                    else:
                        f_bytes = bytes(raw_f)

                    total_len += len(f_bytes)
                    preparsed.add(
                        str(field_name),
                        UploadFile(filename=fname, content=f_bytes, content_type=ctype),
                    )
                if not body_bytes and total_len > 0:
                    body_bytes = b"x" * min(total_len, settings.max_upload_size_bytes + 1)

        req = Request(
            method="POST",
            path=path,
            query_params=q_dict,
            headers=hdrs,
            body=body_bytes,
            preparsed_form=preparsed,
        )
        return ClientResponse(self.app.dispatch(req))


def create_http_handler(app: FastAPI) -> type[BaseHTTPRequestHandler]:
    """Create a standard library `BaseHTTPRequestHandler` bound to `app`."""

    class AppHTTPRequestHandler(BaseHTTPRequestHandler):
        def _handle(self, method: str) -> None:
            parsed = urlsplit(self.path)
            q_dict: Dict[str, str] = {}
            if parsed.query:
                for k, vals in parse_qs(parsed.query, keep_blank_values=True).items():
                    if vals:
                        q_dict[k] = str(vals[-1])

            content_length = int(self.headers.get("Content-Length", "0") or "0")
            body = self.rfile.read(content_length) if content_length > 0 else b""
            hdrs = {k.lower(): v for k, v in self.headers.items()}

            req = Request(
                method=method,
                path=parsed.path or "/",
                query_params=q_dict,
                headers=hdrs,
                body=body,
            )
            resp = app.dispatch(req)

            self.send_response(resp.status_code)
            for hk, hv in resp.finalized_headers().items():
                self.send_header(hk, hv)
            out_bytes = resp.body
            self.send_header("Content-Length", str(len(out_bytes)))
            self.end_headers()
            if method != "HEAD":
                self.wfile.write(out_bytes)

        def do_HEAD(self) -> None:  # noqa: N802
            self._handle("HEAD")

        def do_GET(self) -> None:  # noqa: N802
            self._handle("GET")

        def do_POST(self) -> None:  # noqa: N802
            self._handle("POST")

        def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
            return

    return AppHTTPRequestHandler


def run_http_server(app: FastAPI, host: str = "127.0.0.1", port: int = 8000) -> None:
    """Start a multithreaded HTTP server for `app`."""
    handler_cls = create_http_handler(app)
    server = ThreadingHTTPServer((host, port), handler_cls)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
