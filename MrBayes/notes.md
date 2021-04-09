# Mr. Bayes

First I need to install Mr. Bayes. Claudia had the following in her markdown file:

> 1. Download MrBayes from [here](http://nbisweden.github.io/MrBayes/). In mac:
>
> ```
> brew tap brewsci/bio
> brew install mrbayes --with-open-mpi
> 
> $ which mb
> /usr/local/bin/mb
> ```
>
> Had to troubleshoot a lot!
>
> ```
> brew reinstall mrbayes
> sudo chown -R $(whoami) /usr/local/Cellar/open-mpi/4.1.0
> brew reinstall mrbayes
> ```

I'm going to also try the Brew install command on my own computer.

```sh
brew tap brewsci/bio
brew install mrbayes --with-open-mpi

$ which mb
/usr/local/bin/mb

#unlike claudia I did not have to fix anything
```



