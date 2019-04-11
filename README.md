# rdkit_api

# A Flask API for RDkit

### RDkit, Flask, Swagger for documentation
### Simple authentication with hardcoded username/password for demo purposes only!
### Includes Dockerfile used to deploy on AWS Fargate.

[Demo](https://api.explorablelabs.com) available!

Note that this API uses 2 Dockerfiles: one with RDkit and core libraries pre-installed, and one just for the Flask API (so that changes can be easily made to the Flask API without having to rebuild RDkit). The Dockerfile included here is for the API, and it pulls a private image from hub.docker.com. I have included the Dockerfile for that image in this repo in the dockerhub/ directory.



