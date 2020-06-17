Quick Start 
=============

Source Code
-----------

* `View source code on github <https://github.com/ZekunYin/RabbitSketch>`_

Use RabbitSketch Python APIs
-------------------------------

Install RabbitSketch from pypi:

.. code::

  pip3 install rabbitsketch --user

Install RabbitSketch from source code:

.. code::

  git clone --recursive https://github.com/ZekunYin/RabbitSketch.git
  cd RabbitSketch
  pip3 install . --user

Quick test:

.. code::

  cd RabbitSketch/examples
  python3 test.py genom1.fna genome2.fna

Use RabbitSketch C++ APIs
-------------------------

Install:

.. code::

  git clone --recursive https://github.com/ZekunYin/RabbitSketch.git
  cd RabbitSketch
  mkdir build && cd build
  cmake -DBUILDCXX=ON ..
  make -j4
  make install #todo change it to make install 

Quick test:

.. code::
  
  cd RabbitSketch/examples
  make
  ./test genome1.fna genome2.fna