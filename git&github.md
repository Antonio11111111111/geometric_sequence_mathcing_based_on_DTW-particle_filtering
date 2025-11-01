
# Git 和 GitHub 常用操作指南

## 1. 安装 Git

首先，确保你的电脑已经安装了 Git。访问 [Git 官网](https://git-scm.com/) 下载并安装 Git。

安装完成后，在终端（命令行）输入以下命令检查是否安装成功：

```bash
git --version
````

## 2. 配置 Git

安装 Git 后，需要进行一些基本配置，如设置用户名和邮箱：

```bash
git config --global user.name "你的名字"
git config --global user.email "你的邮箱"
```

## 3. 创建新的 Git 仓库

进入你想管理的项目文件夹，初始化一个新的 Git 仓库：

```bash
git init
```

这会在当前目录创建一个 `.git` 文件夹，表示 Git 开始跟踪该目录下的文件。

## 4. 克隆一个 GitHub 仓库

如果你想要将 GitHub 上的仓库克隆到本地，可以使用以下命令：

```bash
git clone https://github.com/用户名/仓库名.git
```

## 5. 查看仓库状态

查看 Git 仓库的状态，哪些文件被修改了，哪些文件未被追踪等：

```bash
git status
```

## 6. 添加文件到暂存区

将文件添加到 Git 的暂存区，以便提交：

```bash
git add 文件名
```

如果要添加所有更改的文件，可以使用：

```bash
git add .
```

## 7. 提交更改

将暂存区的文件提交到本地仓库：

```bash
git commit -m "提交信息"
```

## 8. 查看提交历史

查看提交的历史记录：

```bash
git log
```

## 9. 创建分支

创建一个新的分支并切换到该分支：

```bash
git checkout -b 分支名
```

## 10. 切换分支

切换到一个已经存在的分支：

```bash
git checkout 分支名
```

## 11. 合并分支

将另一个分支的更改合并到当前分支：

```bash
git merge 分支名
```

## 12. 删除分支

删除一个本地分支：

```bash
git branch -d 分支名
```

## 13. 推送更改到 GitHub

将本地的更改推送到 GitHub 上的远程仓库：

```bash
git push origin 分支名
```

## 14. 拉取远程仓库的更改

拉取 GitHub 上的最新更改到本地仓库：

```bash
git pull origin 分支名
```

## 15. 获取远程仓库的最新提交

如果有其他人提交了更改，可以通过以下命令来更新本地仓库：

```bash
git fetch origin
```

## 16. 查看远程仓库

查看当前远程仓库的配置信息：

```bash
git remote -v
```

## 17. 创建标签（Tag）

在某个特定提交上打标签，通常用于版本发布：

```bash
git tag v1.0
```

## 18. 推送标签到远程仓库

将本地标签推送到 GitHub：

```bash
git push origin v1.0
```

## 19. 删除远程分支

删除远程仓库的分支：

```bash
git push origin --delete 分支名
```

## 20. 回滚到某个历史版本

查看提交记录并找到你需要回滚的提交 ID：

```bash
git log
```

回滚到某个提交：

```bash
git checkout 提交ID
```

如果想彻底回滚并清除修改：

```bash
git reset --hard 提交ID
```

## 21. 解决合并冲突

在合并分支时，如果发生冲突，Git 会标记冲突文件。你需要手动解决冲突，然后重新添加、提交。

查看冲突文件：

```bash
git status
```

打开冲突文件，按照 Git 的标记手动解决冲突。解决完后，将文件标记为已解决：

```bash
git add 冲突文件
```

然后提交：

```bash
git commit -m "解决合并冲突"
```

## 22. GitHub 常用操作

* **Fork 仓库**：在 GitHub 上复制一个仓库到你的 GitHub 账户中。
* **Pull Request**：当你对某个仓库做了更改，可以提交一个 Pull Request 来请求合并更改。
* **Issues**：GitHub 提供的一个功能，用于追踪 bug、任务、功能请求等。

---

> 这些是 Git 和 GitHub 的一些基本操作，希望这份文档能帮助你入门。随着使用的深入，你会逐渐熟悉更多高级命令和操作！


