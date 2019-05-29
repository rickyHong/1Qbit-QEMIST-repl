#!/bin/bash

# Copyright 2019 1QBit
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This is the script that runs on Travis. It simply build the Docker image, runs
# a container, and runs the unit tests inside the container.

docker build $TRAVIS_BUILD_DIR/docker_images/ --tag qemist_test_image

docker run -idt --rm -v $TRAVIS_BUILD_DIR:/root/openqemist --env OMP_THREAD_LIMIT=4 --name qemist_test_container qemist_test_image bash

docker exec -it qemist_test_container bash -c "cd /root/openqemist/openqemist/tests; python3 -m unittest discover"
