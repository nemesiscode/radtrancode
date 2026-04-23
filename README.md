Welcome to the Github repository of the radtrancode project. This project is the development area for the Radtrans and Nemesis codes for radiative transfer calculations and retrievals, originally written in Oxford University.

If you want to use radtrancode/nemesis you have one of two options:
1) You can check out the repository and compile it yourself on your system (details are below)
2) You can pull a pre-compiled 'image' of the code from Docker Hub at https://hub.docker.com/r/patrickirwinoxford/docker_nemesis and run it with Docker. 

If you want to use the docker image, you will also need to add aliases to your shell setup file, which can be found at https://github.com/shubhamvkulkarni/Docker_aliases_for_NEMESIS.git. Depending on whether you use .cshrc, .bashrc or .zsh, you need to download the appropriate add_nemesis_docker_aliases file AND the list of executables, nemesis_executables.txt, into your home directory and run the shell script there. This will add aliases to your shell setup file so that when you run Nemesis, or any other of the radtrancode executables, the system will automatically run the Docker build. All input/output files are passed seamlessly into and out of the Docker 'container'. I am infinitely grateful to my recently graduated D.Phil. student, Shubham Kulkarni, for setting up the Docker route. 
 
If you want to compile radtrancode/nemesis yourself (and thus be able to edit the code and recompile as you wish), please do the following steps: 

-------------------------

On your system go to the directory where you want to put the radtrancode git repository, i.e.,

cd somewhere safe, for example

```
cd ~/repositories
```

then type

```
git clone https://github.com/nemesiscode/radtrancode.git
```

You should now have the current git repository in your safe directory.

Now on your system type the command:

```
git config --global --list
```

If you don’t see something like:

```
user.name=Firstname Lastname						{obviously your name here}
user.email=firstname.lastname@...             {obviously your email address here}
```

then you should add them before checking anything in or out. To do this type:

```
git config --global user.name “Firstname Lastname" 
git config --global user.email “firstname.lastname@...  “ 
```

Also, if you want to change the default git editor from “vi” or “vim” type something like:
`git config --global core.editor “pico"`

Finally, if you want the system to give you colour-coded git status messages (which do look quite pretty) type:

```
git config --global color.ui "auto"
```

To update your copy with the latest version on github the command is:

```
git pull origin
```

(updates the default branch master)

To check in your updates to github you type

```
git push origin <branchname>
```
To compile the code, please look at the instrunctions in AACOMPILE.txt. Remember to check the ``ISYS`` variable for your specific compilation (check AACOMPILE.txt).
It is also helpful to add the following to your BASHRC/CSHRC file:

``
export RADREPO=/PATH-TO/radtrancode/
``

...as this will prevent you having to update the `dataarchive.f` code every time to point
to your `raddata` directory.
