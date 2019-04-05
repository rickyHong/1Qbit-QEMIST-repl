#!/bin/bash

# This is the script that runs on Travis. It simply build the Docker image, runs
# a container, and runs the unit tests inside the container.

docker build $TRAVIS_BUILD_DIR/docker_images/ --tag qemist_test_image

docker run -idt --rm -v $TRAVIS_BUILD_DIR:/root/QEMIST --env OMP_THREAD_LIMIT=4 --name qemist_test_container qemist_test_image bash

docker exec -it qemist_test_container bash -c "cd /root/QEMIST/qemist/tests; python3 -m unittest discover"
