Contributing
=============

Contributions, either with `issues <https://github.com/pierre-24/pymerk/issues>`_ or `pull requests <https://github.com/pierre-24/pymerk/pulls>`_ are welcomed.

Install
-------

If you want to contribute, this is the usual deal: start by `forking <https://guides.github.com/activities/forking/>`_, then clone your fork and use 
the following install procedures instead.

.. code-block:: bash

    cd pymerk

    # definitely recommended in this case: use a virtualenv!
    python -m venv virtualenv
    source virtualenv/bin/activate

    # install also dev dependencies
    make install-dev

Tips to Contribute
------------------

* **Start with an issue**: A good place to start is the `list of issues <https://github.com/pierre-24/pymerk/issues>`_.
  In fact, it is easier if you start by filling an issue, and if you want to work on it, say so there, so that 
  everyone knows that the issue is being handled.

* **Work on a separate branch**: Don't work directly on ``main``.
  Since this project follows the `git flow branching model <http://nvie.com/posts/a-successful-git-branching-model/>`_,
  you should base your branch on ``main``:

  .. code-block:: bash

      git checkout -b my_feature origin/main

* **Run linting and tests regularly**:

  .. code-block:: bash

      make lint
      make test

  The code follows the `PEP-8 style recommendations <http://legacy.python.org/dev/peps/pep-0008/>`_, 
  checked by `flake8 <https://flake8.pycqa.org/en/latest/>`_.
  Having an extensive test suite is also a good idea to prevent regressions.

* **Pull requests should be**:

  - **Unitary**: Focus on a single feature or fix per PR
  - **Tested**: Include unit test(s) for the change
  - **Documented**: Update docstrings and documentation if needed
  - **Passing**: Both the test suite and linting must pass

  All checks must succeed for the pull request to be accepted.