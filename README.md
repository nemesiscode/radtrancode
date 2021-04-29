Welcome to the new Github repository of the radtrancode project. This project is the development area for the Radtrans and Nemesis codes for radiative transfer calculations and retrievals, originally written in Oxford University.



If you are new to the system and want to check out the radtrancode repository onto you local system, please do the following steps: 

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
git config --global core.editor “pico"

Finally, if you want the system to give you colour-coded git status messages (which do look quite pretty) type:

```
git config --global color.ui "auto"
```

To update your copy with the latest version on github the command is:

```
git pull origin
```

(updates the default branch master)

To check in your updates to github/gitlab you type

```
git push origin <branchname>
```
