Linux Primer
============

Linux path environment variable
-------------------------------

(adapted from https://linuxconfig.org/linux-path-environment-variable)

Using Linux ``PATH`` variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Linux path environmental variable contains a list of directories in which the shell looks for executable programs every time you run a command or program.
Use ``echo`` command to print your ``PATH`` variable:

::

    $ echo $PATH
    /home/lilo/bin:/usr/local/bin:/usr/bin:/bin:/usr/games

If the program / command is located within the ``PATH``, user do not need to include full path in order to execute a certain command.
For example: ``date`` command is located within ``/bin``:

::

    $ which date
    /bin/date

and ``/bin`` is defined in the ``PATH`` variable.
Therefore, to execute date command is easy as:

::

    $ date

Adding a new directory into ``PATH`` variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From time to time you may need to add new directory into your ``PATH`` scope.
The following example adds new directory ``/bin/myscripts`` into ``PATH`` variable:

::

    $ echo $PATH
    /home/lilo/bin:/usr/local/bin:/usr/bin:/bin:/usr/games
    $ PATH=$PATH:/bin/myscripts
    $ export PATH
    $ echo $PATH
    /home/lilo/bin:/usr/local/bin:/usr/bin:/bin:/usr/games:/bin/myscripts

The above method temporarily adds certain directories into the ``PATH`` variable.
To permanently modify the ``PATH`` variable, copy following text to ``/home/yourname/.bash_profile`` or ``/home/yourname/.bashrc``:

::

    PATH=$PATH:/bin/myscripts
    export PATH

And then execute this command: ``source /home/yourname/.bash_profile`` or ``source /home/yourname/.bashrc``, depending on which file you have modified.

To see whether the ``PATH`` has been successfully modified, execute: ``echo $PATH``.
