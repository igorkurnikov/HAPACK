FROM nsimakov/harlem_ready:1


LABEL description="image to run tests with git repo location like in shippable"

# copy repo
VOLUME /root/src/gitlab.com/mkurnikovagroup/HAPACK

