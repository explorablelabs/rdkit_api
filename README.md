# rdkit_api

# A Flask API for RDkit

* RDkit, Flask, Swagger for documentation
* Simple authentication with hardcoded username/password for demo purposes only!
* Includes Dockerfile used to deploy on AWS Fargate.

### [Demo](https://api.explorablelabs.com) here!

* Docker: Note that this API uses 2 Docker images: one with RDkit and core libraries pre-installed, and one just for the Flask API (so that changes can be quickly made to the API without having to rebuild RDkit). The Dockerfile included at the top level here is for the API only, and it pulls a private image from hub.docker.com. You can find the Dockerfile for that image in this repo in the dockerhub/ directory. You can probably use this API with just the image informaticsmatters/rdkit-python-debian. The reason for creating a separate image here was to be able to deploy with CircleCI.



