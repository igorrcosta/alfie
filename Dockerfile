FROM python:3.10
ADD *.py ./
ADD test/ ./test/
ADD bin/* /usr/local/bin/
ADD out/ ./out/
ADD in/ ./in/
RUN pip install --no-cache-dir biopython numpy
CMD ["python", "alfie.py", "-r", "test/homo_y.fasta", "-i", "test/pan_y.fasta", "-g" ,"test/y.gtf"]

