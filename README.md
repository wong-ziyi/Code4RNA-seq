# Code4RNA-seq
Code book for RNA-seq analysis (From FASTQ to Statistical Results)  
Check [here](https://github.com/wong-ziyi/Code4RNA-seq/blob/master/00.InitialSetUp4WSLonWindows.md) to set up a Windows PC.
## Use git on Ubuntu terminal
### Installation
```bash
$ sudo apt-get install git
```
### Clone this Git repository into a specific folder and add remote
```bash
$ sudo git clone https://github.com/wong-ziyi/Code4RNA-seq /mnt/g/Code4RNA-seq
$ cd /mnt/g/Code4RNA-seq
$ git remote add origin https://github.com/wong-ziyi/Code4RNA-seq
```
### Update your local repository
```bash
$ sudo git pull
```
### Push a local commit
```bash
$ touch test.txt
$ git add test.txt
$ git commit -m "test"
$ git push
```
### Delete a file
```bash
$ git rm test.txt
$ git commit -m "test2"
$ git push
```
### Please see more detials on [here](https://dont-be-afraid-to-commit.readthedocs.io/en/latest/git/commandlinegit.html)
