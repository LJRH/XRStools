Installation
============

  * If you install from a Debian package you can skip the following points, install it , and then go directly to the code invocation section

  * Using Git, sources can be retrived with the following command ::

      git clone https://gitlab.esrf.fr/mirone/XRStools

  * for a local installation you can use  ::

	python setup.py install --prefix=~/packages

    then to run the code you must do beforehand ::

	export PYTHONPATH=/home/yourname/packages/lib/python2.7/site-packages
	export PATH=/home/yourname/bin:$PATH


  * For the roi selection tool you need a recent version of pymca installed on your sistem


  * Examples :
        There is an archive containing datas for tests by anonymous ftp at ::

	  http://ftp.esrf.fr/pub/scisoft/___STILL_TO_BE_PLACED
	  
