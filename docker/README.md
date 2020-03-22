# Ubuntu based Docker image for Harlem testing

Here are the notes on the preparation of docker container with dependencies for Harlem 
installation and testing.

## Creating Image

### Making Image

```shell
REM make image
cd C:\MYPROG\HAPACK\docker
docker build -t nsimakov/harlem_ready:2 -f Dockerfile .
```

### Run

```bash
#run
docker run --name harlemtest -h harlemtest --rm -it nsimakov/harlem_ready:1

### Push

```bash
#push to docker cloud
docker push nsimakov/harlem_ready:2
```

## Testing Harlem In Docker Container Locally

```shell
REM run
docker run -v /C/MYPROG/HAPACK:/root/src/gitlab.com/mkurnikovagroup/HAPACK ^
          --name harlemtest -h harlemtest -it --rm nsimakov/harlem_ready:2 bash
```

or without setting names

```shell
REM run
docker run -v /C/MYPROG/HAPACK:/root/src/gitlab.com/mkurnikovagroup/HAPACK -it --rm nsimakov/harlem_ready:2 bash
```
