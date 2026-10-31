# Getting started with dndsR

This tutorial installs dndsR on Linux with Singularity and confirms that the command line works. It uses the same tested commands as the AI agent skill. It does not run an analysis yet, because validated example data is still being prepared. See the [README](../../README.md) for the documented analysis steps.

You need bash, git, Singularity, and network access. The image is 1.3 GB.

## 1. Get the code

```
git clone https://github.com/Scrumpis/dndsR
cd dndsR
```

## 2. Get the container

This exact digest always gives the same image.

```
singularity pull dndsr.sif docker://scrumpis/dndsr@sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831
```

## 3. Install the `dndsr` command

This installs into your home directory (`~/.config/dndsr`, `~/.dndsr/rlib`, and `~/.local/bin`), adds shell completions, and may edit `~/.zshrc`. Keep the cloned folder.

```
tools/dndsr-install install --engine singularity --sif dndsr.sif
```

If it says `~/.local/bin` is not on your `PATH`, run the line it prints.

## 4. Check that it works

```
tools/dndsr-install doctor
dndsr --help
```

You should see 11 commands, including `calculate_dnds`. For help with a command:

```
dndsr calculate_dnds --help
```
