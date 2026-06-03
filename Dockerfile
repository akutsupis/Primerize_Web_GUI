# 1. Use an optimized, official Python lightweight image
FROM python:3.11-slim

# 2. Prevent Python from writing .pyc files and buffering stdout/stderr
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /code

# 3. Copy only requirements first to leverage Docker caching layers
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# 4. Copy the rest of your application code (including the primerize folder)
COPY . .

# 5. Hugging Face Spaces require Docker containers to listen on port 7860
EXPOSE 7860

# 6. Run Streamlit and bind it to the specific port and address routing HF expects
CMD ["streamlit", "run", "app.py", "--server.port=7860", "--server.address=0.0.0.0"]