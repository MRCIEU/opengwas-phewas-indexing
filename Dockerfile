FROM python:3.11-alpine

COPY ./requirements.txt /requirements.txt
RUN python -m pip install -r /requirements.txt

COPY ./phewas-indexing /phewas-indexing

ENTRYPOINT ["python", "/phewas-indexing/main.py"]
