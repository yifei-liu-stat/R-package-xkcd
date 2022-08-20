# xkcd
Developed an R package and tested it for STAT8054 [Doing Assignment 2](http://www.stat.umn.edu/geyer/8054/hw/).

This is a test R package that is used for STAT8054 [Doing Assignment 2](http://www.stat.umn.edu/geyer/8054/hw/). A rough idea about xkcd distribution is shown in [xkcd comic 2118](https://xkcd.com/2118/). More details can be found in the assignment description section.

To build and check this package, run the following in the parent directory containing `xkcd` in `bash`:

```bash
R CMD build xkcd
R CMD check xkcd_0.0.0.9000.tar.gz
```

This package also supports C implementation. This can be controlled by the `flag` argument (default is `flag = "R"`) in each core function (`rxkcd`, `dxkcd`, `pxkcd`, `qxkcd`). For example, 

```R
rxkcd(1000, nsd = 2, flag = "C")
dxkcd(0.1, swap.end.points = TRUE, flag = "C")
pxkcd(0.1, log.p = TRUE, flag = "C")
qxkcd(0.5, flag = "C")
```




