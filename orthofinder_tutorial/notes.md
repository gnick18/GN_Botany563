# Orthofinder tutorial Notes

## Installation and running on the example dataset

I first installed the program using wget and bioconda using the following.

```shell
wget https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder.tar.gz
tar xzvf OrthoFinder.tar.gz 
conda install -c bioconda orthofinder
```

The installation step with Bioconda took a very long time. However when i downloaded the data set tar link it also installed the full program.

To check that it was installed correctly one can use `orthofinder -h`

this worked!

I will now run the analysis on the example data set

```shell
orthofinder -f ExampleData/
```

This put the results inside the ExampleData folder under `/OrthoFinder/Results_Feb23/`



