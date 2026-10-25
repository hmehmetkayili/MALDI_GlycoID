# Deployment & Security Guide — MALDI-GlycoID

## 1. Local Deployment
```bash
python scripts/run_local.py --host 127.0.0.1 --port 8000
```

## 2. Online Deployment on Render.com (via GitHub)
The repository includes pre-configured Render deployment manifests (`render.yaml`, `Procfile`, and `runtime.txt`) and automatic `PORT` / `0.0.0.0` binding in `app/config.py`.

### Option A: One-Click Blueprint Deployment (`render.yaml`)
1. Push this repository to your GitHub account.
2. Log in to [Render.com](https://dashboard.render.com/) and click **New +** &rarr; **Blueprint**.
3. Connect your GitHub repository (`MALDI-GLYCOID2`).
4. Render automatically detects `render.yaml`, installs `requirements.txt`, starts `python scripts/run_local.py --host 0.0.0.0`, and verifies `/api/health`.

### Option B: Manual Web Service Setup on Render
1. On Render Dashboard, click **New +** &rarr; **Web Service** and select your GitHub repository.
2. Configure the service settings:
   - **Runtime:** `Python 3` (or `Docker`)
   - **Build Command:** `pip install -r requirements.txt`
   - **Start Command:** `python scripts/run_local.py --host 0.0.0.0`
3. Click **Create Web Service**. Render automatically injects the `PORT` environment variable and brings your URL online.

## 3. Containerized Docker Deployment
```bash
# Build and start via Docker Compose
docker compose up --build -d

# Check container health
docker ps
```
The container runs as a non-root user (`glycoid`, UID 1000) and includes an automated `HEALTHCHECK` executing `python scripts/run_local.py --check-only`.

## 3. Built-In Security Controls
1. **Strict Upload Size Limit:** Enforces `MALDI_GLYCOID_MAX_UPLOAD_BYTES` (default 50 MB) with HTTP `413` rejection for oversized payloads.
2. **Filename Sanitization & Path Traversal Protection:** Rejects `..`, leading slashes, and drive prefixes in both single-file uploads (`/api/spectra/upload`) and batch `.zip` archives (`extract_batch_zip_safely`).
3. **ZIP Bomb & Nested Archive Rejection:** Rejects nested archives (`.zip`, `.tar`, `.gz`, `.7z`, `.rar`), enforces maximum file count (`MALDI_GLYCOID_MAX_BATCH_FILES`), and limits total uncompressed byte size in memory.
4. **HTTP Security Headers:** Every response includes `X-Content-Type-Options: nosniff`, `X-Frame-Options: DENY`, `Referrer-Policy: strict-origin-when-cross-origin`, `X-XSS-Protection: 1; mode=block`, and `Content-Security-Policy`.
