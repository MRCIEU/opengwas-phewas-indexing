FROM python:3.11-alpine

COPY ./phewas-indexing/requirements.txt /phewas-indexing/requirements.txt
RUN python -n pip install -r /phewas-indexing/requirements.txt

COPY ./phewas-indexing /phewas-indexing

ENTRYPOINT ["python", "/phewas-indexing/main.py"]
