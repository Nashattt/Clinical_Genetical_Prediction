# Use the official Python 3.12.3 slim image
FROM python:3.12.3-slim

# Set the working directory in the container
WORKDIR /app

# Copy the dependencies file
COPY requirements.txt .

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application code
COPY . .

# Expose port 7860 for Hugging Face Spaces
EXPOSE 7860

# Run with Gunicorn
CMD ["gunicorn", "--bind", "0.0.0.0:7860", "app:app"]