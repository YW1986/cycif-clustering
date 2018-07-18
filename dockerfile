FROM python:3.6
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt && \
	mkdir /input && \
	mkdir /output
VOLUME /input
VOLUME /output
VOLUME /config
COPY cycif.py /usr/local/bin/cycif.py
ENTRYPOINT ["python", "/usr/local/bin/cycif.py"]