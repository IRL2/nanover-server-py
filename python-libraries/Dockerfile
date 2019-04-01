FROM continuumio/anaconda3
ADD narupa-protocol /
RUN pip install -r requirements.txt
ENV PYTHONPATH "${PYTHONPATH}:/:/narupa:/narupy"
WORKDIR /narupy/benchmark
EXPOSE 8000
EXPOSE 8001
CMD ["python", "benchmark_server.py"]