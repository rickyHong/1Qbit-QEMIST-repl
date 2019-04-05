# Examples

If running Jupyter notebooks inside the docker container, use run the docker container with the port that Jupyter is serving on exposed. The default is 8888, so use the `-p 8888:8888` argument.
Then run the Jypiter notebook with
```bash
$jupyter notebook --no-browser --allow-root --ip=0.0.0.0
```
