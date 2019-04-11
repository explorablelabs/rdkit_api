FROM kathrynloving/explorable-labs:0.0.1

COPY . /app
EXPOSE 5000

WORKDIR /app
RUN pip install -r requirements.txt
ENV PYTHONPATH "${RBASE}:${PYTHONPATH}"
ENTRYPOINT ["python"]
CMD ["app.py"]
