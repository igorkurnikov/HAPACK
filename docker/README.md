# Centos-7 based Docker image for Harlem testing

Here are the notes on the preporation of docker container with dependencies for Harlem installation and testing.

## Creating Image

### Making Slurm RPMs

First we need slurm RPMs.
DockerfileMakeSlurmRPM describes simple image for centos 7 rpm making.
Here is listing on the whole process

```bash
# make image
docker build -t nsimakov/harlem_ready:1 -f Dockerfile .
```

### Run

```bash
#run
docker run --name harlemtest -h harlemtest --rm -it nsimakov/harlem_ready:1

### Push

```bash
#push to docker cloud
docker push nsimakov/harlem_ready:1
```

## Testing Harlem In Docker Container Locally

```bash
#make image in akrr root
docker build -t pseudo_repo/harlemtest:latest .

#run
docker run -v /C/MYPROG/HAPACK:/root/src/gitlab.com/mkurnikovagroup/HAPACK \
          --name harlemtest -h harlemtest -it --rm pseudo_repo/harlemtest:latest bash
```