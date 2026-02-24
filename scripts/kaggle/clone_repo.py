"""Git repository cloning utilities for Kaggle notebooks.

Clone repositories with optional authentication and branch selection.
"""

import subprocess
import os
from pathlib import Path


def clone_repo(
    repo_url,
    dest_path="/kaggle/working/repo",
    branch=None,
    token=None,
    username=None
):
    """Clone a git repository to a local path.

    Args:
        repo_url (str): Repository URL (HTTP or SSH)
        dest_path (str): Local destination path. Default /kaggle/working/repo
        branch (str, optional): Specific branch to checkout
        token (str, optional): GitHub personal access token for private repos
        username (str, optional): GitHub username (used with token)

    Returns:
        Path: Path to cloned repository

    Example:
        >>> # Public repo
        >>> repo = clone_repo("https://github.com/user/repo.git")

        >>> # Private repo with token
        >>> repo = clone_repo(
        ...     "https://github.com/user/private-repo.git",
        ...     token="ghp_xxxx",
        ...     username="user"
        ... )

        >>> # Specific branch
        >>> repo = clone_repo(
        ...     "https://github.com/user/repo.git",
        ...     branch="development"
        ... )
    """
    dest = Path(dest_path)

    # If destination exists and is a git repo, skip cloning
    if dest.exists() and (dest / ".git").exists():
        print(f"  ⚠ Repository already exists at {dest}")
        if branch:
            print(f"  → Checking out branch: {branch}")
            subprocess.run(
                ["git", "-C", str(dest), "checkout", branch],
                check=True
            )
        return dest

    print(f"=== Cloning repository ===")
    print(f"  URL: {repo_url}")
    print(f"  Dest: {dest}")

    # Build authenticated URL if token provided
    clone_url = repo_url
    if token and username:
        # Convert https://github.com/user/repo.git to https://user:token@github.com/user/repo.git
        if repo_url.startswith("https://github.com/"):
            clone_url = repo_url.replace(
                "https://github.com/",
                f"https://{username}:{token}@github.com/"
            )

    # Clone command
    cmd = ["git", "clone"]
    if branch:
        cmd.extend(["-b", branch])
    cmd.extend([clone_url, str(dest)])

    subprocess.run(cmd, check=True)
    print(f"  ✓ Cloned successfully")

    return dest


def pull_latest(repo_path, branch=None):
    """Pull latest changes from remote.

    Args:
        repo_path (str): Path to local repository
        branch (str, optional): Branch to pull. If None, pulls current branch.

    Example:
        >>> pull_latest("/kaggle/working/repo", branch="main")
    """
    print(f"=== Pulling latest changes ===")

    if branch:
        subprocess.run(
            ["git", "-C", repo_path, "checkout", branch],
            check=True
        )

    subprocess.run(
        ["git", "-C", repo_path, "pull"],
        check=True
    )
    print(f"  ✓ Updated to latest")


def get_repo_info(repo_path):
    """Get information about a cloned repository.

    Args:
        repo_path (str): Path to local repository

    Returns:
        dict: Repository info (branch, commit, remote)

    Example:
        >>> info = get_repo_info("/kaggle/working/repo")
        >>> print(f"Current branch: {info['branch']}")
    """
    info = {}

    # Get current branch
    result = subprocess.run(
        ["git", "-C", repo_path, "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True, text=True, check=True
    )
    info["branch"] = result.stdout.strip()

    # Get current commit
    result = subprocess.run(
        ["git", "-C", repo_path, "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True, check=True
    )
    info["commit"] = result.stdout.strip()

    # Get remote URL
    result = subprocess.run(
        ["git", "-C", repo_path, "remote", "get-url", "origin"],
        capture_output=True, text=True, check=True
    )
    info["remote"] = result.stdout.strip()

    return info


def setup_repo(
    repo_url,
    dest_path="/kaggle/working/docking-pipeline",
    branch=None,
    token=None,
    username=None,
    pull_if_exists=False
):
    """Complete repository setup: clone or update.

    Args:
        repo_url (str): Repository URL
        dest_path (str): Local destination
        branch (str, optional): Specific branch
        token (str, optional): GitHub token for private repos
        username (str, optional): GitHub username
        pull_if_exists (bool): Pull latest if repo exists. Default False.

    Returns:
        dict: Repository info and path

    Example:
        >>> repo = setup_repo(
        ...     "https://github.com/user/docking.git",
        ...     dest_path="/kaggle/working/docking-pipeline",
        ...     branch="main"
        ... )
        >>> print(repo['path'])
        >>> print(repo['info']['commit'])
    """
    dest = Path(dest_path)

    # Clone or update
    if dest.exists() and (dest / ".git").exists():
        if pull_if_exists:
            pull_latest(str(dest), branch)
        else:
            print(f"  ✓ Repository already exists at {dest}")
    else:
        clone_repo(repo_url, str(dest), branch, token, username)

    # Get info
    info = get_repo_info(str(dest))

    print(f"\nRepository ready:")
    print(f"  Path: {dest}")
    print(f"  Branch: {info['branch']}")
    print(f"  Commit: {info['commit']}")

    return {"path": dest, "info": info}
