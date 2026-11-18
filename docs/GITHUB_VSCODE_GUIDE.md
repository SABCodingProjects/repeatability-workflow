# Upload this project to GitHub from VS Code

You already have VS Code and a GitHub account. These steps assume you want to create a **new repository** for this project.

## 1. Open the extracted project

1. Extract the ZIP file.
2. In VS Code, select **File → Open Folder**.
3. Select the `repeatability-workflow` folder containing `README.md`, `config.R` and `run_analysis.R`.

Open this folder itself, not its parent folder or your whole research directory.

## 2. Check that Git is available

Select **Terminal → New Terminal**, then run:

```sh
git --version
```

If a version appears, continue. If Git is not found, install [Git](https://git-scm.com/downloads) and restart VS Code. Having a GitHub account does not necessarily mean Git is installed locally.

## 3. Initialise the repository

1. Click **Source Control** in the left sidebar (the branching icon).
2. Select **Initialize Repository**.
3. Check that it initialises the `repeatability-workflow` folder.

If this folder is already a repository, do not initialise another one; continue to reviewing its files.

## 4. Review and commit the files

1. Under **Changes**, click the files to review them.
2. Confirm that only generic code, documentation, tests and synthetic examples are present.
3. Real data belong in `data/private/`, and analysis outputs belong in `outputs/`. Both are ignored by Git. Keep the public `config.R` set to the synthetic example.
4. Click **+** beside each intended file, or **Stage All Changes** after checking the entire list.
5. Enter a commit message, for example: `Add commented recovery repeatability workflow`.
6. Click **Commit**.

If Git asks for your identity, run these commands in the terminal, replacing the placeholders, then retry the commit:

```sh
git config user.name "Your chosen name"
git config user.email "Your GitHub email or provided noreply address"
```

GitHub's **Settings → Emails** shows your available noreply address. These local commands set the author recorded in this repository's commits; they do not sign you into GitHub.

## 5. Publish to GitHub

1. Open the Command Palette: **Ctrl+Shift+P** on Windows/Linux or **Cmd+Shift+P** on macOS.
2. Search for **Publish to GitHub** and select it.
3. Sign in with your existing account if prompted.
4. Enter a repository name, such as `recovery-repeatability-workflow`.
5. Choose **private** while reviewing, or **public** when ready to share.
6. Review any file-selection prompt and complete publication.
7. Select **Open on GitHub** to view the repository.

The command creates the GitHub repository and uploads your committed files. You do not need to create an empty repository on the website first.

## 6. Verify the upload

On GitHub:

1. Check that the README displays correctly and only synthetic data are present.
2. Open **Actions → R workflow checks** and check whether the tests pass.
3. If tests fail, read the log, correct the problem locally, and commit and push the fix.

The scripts were not executed during preparation because R was unavailable. Run the checks locally as well if you have R installed:

```sh
Rscript scripts/install_dependencies.R
Rscript tests/run_tests.R
```

Passing checks verifies the synthetic workflow; it does not establish that a model is suitable for every research dataset.

## Later updates

Edit → review the changed files → stage → commit → **Push** or **Sync Changes** in Source Control. Check the Actions result again.

Do not upload research data just to demonstrate the code. Ignore rules do not remove files already committed to Git history.

Official reference: [VS Code source control quickstart](https://code.visualstudio.com/docs/sourcecontrol/quickstart).
