FROM python:3.11-slim

LABEL maintainer="MALDI-GlycoID Scientific Software Team"
LABEL description="Scientifically transparent MALDI-TOF MS N-glycan profiling and glycopeptide screening platform"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLCONFIGDIR=/tmp/maldi_glycoid_mpl \
    MALDI_GLYCOID_HOST=0.0.0.0 \
    MALDI_GLYCOID_PORT=8000

WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY . /app

# Create non-root user for secure container execution
RUN useradd -m -u 1000 glycoid && chown -R glycoid:glycoid /app
USER glycoid

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python scripts/run_local.py --check-only || exit 1

CMD ["python", "scripts/run_local.py", "--host", "0.0.0.0"]
