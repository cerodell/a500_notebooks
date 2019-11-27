# How to fork a git repo


These instructions are a combination of individual work by:

Chase Pettit
https://gist.github.com/Chaser324/ce0505fbed06b947d962

&

Phil Austin
https://github.com/phaustin/a500_notebooks/blob/master/notebooks/course_helpfiles/Readme_github.md

##### Note that a500_notebooks (by Phil Austin) is the repo thats being forked/cloned


## Creating a Fork
Go to GitHub login to you account.

https://github.com

Go to Phil Austins a500_notebooks repository and click the `fork` button. 

https://github.com/phaustin/a500_notebooks


<!-- #region -->
### Clone the fork to your local machine

###### (HTTPS)
`git clone https://github.com/<$USER>/a500_notebooks.git`

###### or


###### (ssh)
`git clone git@github.com:<$USER>/a500_notebooks.git`



##################################################################

<!-- #endregion -->

<!-- #region -->
## Add upstream repo to a Fork

######  Change into cloned repo directory

`cd a500_notebooks/`

### Add 'upstream' repo to list of remotes

###### (HTTPS)
`git remote add upstream https://github.com/phaustin/a500_notebooks.git`


###### or


###### (ssh)
`git remote add upstream git@github.com:phaustin/a500_notebooks.git`


### Verify the new remote named 'upstream'

`git remote -v`


##################################################################

<!-- #endregion -->

<!-- #region -->
## Working with a Fork



Whenever you want to update your fork with the latest upstream changes, you'll need to first fetch the upstream repo's branches and latest commits to bring them into your repository:

### Fetch from upstream remote
`git fetch upstream`


##################################################################

<!-- #endregion -->

<!-- #region -->
## Creating a branch within a Fork


### Create a new branch with your initials
`git checkout -b <initials>`

### View all branches, including those from upstream
`git branch -va`

##################################################################

<!-- #endregion -->

<!-- #region -->
## Working with a branch

### Checkout your master branch and merge upstream
`git checkout master`

`git merge upstream/master`


### Push your branch up to your fork

`git push origin <initials>`

### Issue a pull request by checking out your branch on github and hitting the pull request button

##################################################################

##################################################################

<!-- #endregion -->

# Day to day workflow

<!-- #region -->
## When changes are made on the upstream master, you can incorporate those changes to your fork like this:

Fetch project branches from the upstream repository to get all the commits. Your commits to master will be stored in the local branch upstream/master.

### Fetch new upstream changes


` git fetch upstream `

###  Checkout the master branch from your local fork.


` git checkout master `

### Now merge the changes from upstream/master into your local master branch.

Your fork’s master branch will be in sync with the upstream repository. You will not lose your local changes.

NOTE: You will be prompted to type a commit

`git merge upstream/master`

##################################################################

<!-- #endregion -->

<!-- #region -->
## Once you have incorporated changes to you local master you can merge your local branch


### Make sure you are still in the master branch 

` git checkout master `

### Add commit on new changes

`git commit -m "Add note" `


### Merge your local branch to the local master 

NOTE: You will be prompted to type a commit

`git merge <initials>`

### Push you local master with you merged local branch 

`git push`



###### Also useful for updating your local master from your branch (see below for update to branch)


##################################################################
<!-- #endregion -->

<!-- #region -->
## Workflow for your local branch repo:

### To add new changes
`git add .`


### Add commit on new changes

`git commit -m "Add Note"` 


### Push changes to branch

` git push origin <initials> `
<!-- #endregion -->

```python

```
