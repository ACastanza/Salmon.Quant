### copyright 2017-2021 Regents of the University of California and the Broad Institute. All rights reserved.
FROM genepattern/salmon-quant\:beta

MAINTAINER Barbara Hill <bhill@broadinstitute.org>

# docker build --rm https://github.com/genepattern/Salmon.Quant.git#develop -f Dockerfile -t genepattern/salmon-quant:<tag>
# make sure this repo and tag match the manifest & don't forget to docker push!
# docker push genepattern/salmon-quant:<tag>

# you can use this command to run Docker and iterate locally (update for your paths and module name, of course)
# docker run --rm -it -v /c/Users/MyUSER/PathTo/Salmon.Quant:/mnt/mydata:rw genepattern/Salmon.Quant:<tag> bash