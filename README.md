# Visflux

This a local python package based on d3flux, aimed for visualization of cobra models, providing basic interactive functionalities.

![](https://vkceyugu.cdn.bspapp.com/VKCEYUGU-imgbed/b53d237f-4de9-49ef-88f5-78904d7f1e01.png)

*Note*: This is still very much alpha software. Many things won't work :)

## Install

1. Download this repo as a zip file. Unpack it and navigate to the root directory.
2. Install via pip locally:

```
pip install .
```

3. Go to the directory ./visflux/main, copy file "create_map.py" to a desired directory with a COBRA model in JSON format, then execute the python script. Example:

```
ls
$ create_map.py model.json
python ./create_map.py ./model.json
```

Now you will get a single html file named "model.html" in the same directory. Open it in a browser and you shall see the visualization result. Note that an Internet connection is required since this webpage uses CDN services.

## Documentation for D3flux

[https://pstjohn.github.io/d3flux/](https://pstjohn.github.io/d3flux/)