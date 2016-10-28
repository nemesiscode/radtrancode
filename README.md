Welcome to the new git repository of the radtrancode project. This project is the development area for the Radtrans and Nemesis codes for radiative transfer calculations and retrievals, written in Oxford.



If you are new to the system and want to check out the radtrancode repository onto you local system, please do the following steps: 

-------------------------

First create a ‘key pair’ on your system. To do this (Linux or OSX) type:

ssh-keygen –t rsa

It will default to writing to ~/.ssh/id_rsa. If a key already exists it will ask if you want to overwrite this or make a new key pair somewhere else.  It will then ask you for a passphrase. Please enter one to make the system more secure. The public part of the key is written (default) to ~/.ssh/id_rsa.pub

Now log on to the physics gitlab system at
https://gitlab.physics.ox.ac.uk

Click on the icon at top right to access ‘profile settings’. Go to ‘SSH keys’. Copy and paste the contents of to ~/.ssh/id_rsa.pub into the box, give the key a title such as ‘Key on oxpln20’ and then click ‘add key’. 

Now, on your system go to the directory where you want to put the radtrancode git repository, i.e.

cd somewhere safe

then type

git clone git@gitlab.physics.ox.ac.uk:planetary/radtrancode.git

You should now have the current git repository in your safe directory.

Now on your system type the command:
git config --global --list

If you don’t see something like:

user.name=Patrick Irwin						{obviously your name here}

user.email=patrick.irwin@physics.ox.ac.uk             {obviously your email address here}

then you should add them before checking anything in or out. To do this type:

git config --global user.name “Patrick Irwin" 
git config --global user.email “patrick.irwin@physics.ox.ac.uk “ 


Also, if you want to change the default git editor from “vi” or “vim” type something like:
git config --global core.editor “pico"

Finally, if you want the system to give you colour-coded git status messages (which do look quite pretty) type:

git config --global color.ui "auto"

To update your copy with the latest version on github the command is:

git pull origin
(updates the default branch master)

To check in your updates to github you type
git push origin <branchname>

If you are a git 'master' you can push your copy of the master branch onto the repository. If you are a 'developer' you can only push branches with a different name. I can then choose to merge them with the master branch from within github.

Note that if you get weird arcane error messages when you try to 'push' it's probably because the repository has since changed and you need to type 'git pull origin' first to make sure everything is in sync.
