FROM python:3.6
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir Cython && \
	pip install --no-cache-dir -r /tmp/requirements.txt && \
	pip uninstall -y Cython
ADD cycif.py /
ENTRYPOINT ["python", "./cycif.py"]
