# Choose an appropriate base image
FROM hannahdoerpholz/radarplot:latest

RUN useradd -ms /bin/bash newuser


# Set the working directory inside the container
WORKDIR /app

# Copy the Conda packages list to the container
#COPY requirements.txt .
COPY radar_charts.py /app/
COPY inputs/ /app/inputs/
COPY outputs/ /app/outputs/

RUN chown -R newuser:newuser /app/inputs/
RUN chown -R newuser:newuser /app/outputs/
RUN chown newuser:newuser /app/radar_charts.py

USER newuser
# Entrypoint command 
CMD ["python", "radar_charts.py"]
