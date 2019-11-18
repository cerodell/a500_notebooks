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

### Fetch new upstream changes


` git fetch upstream `

###  Switch (checkout) your branch

` git checkout -b <initials> `


` git rebase upstream/master `

##################################################################

<!-- #endregion -->

<!-- #region -->
## Workflow for your forked/branch repo:

### To add new changes
`git add .`


### Add commit on new changes

`git commit -m "Add Note"` 


### Push changes to branch

` git push origin <initials> `
<!-- #endregion -->

```python

```
