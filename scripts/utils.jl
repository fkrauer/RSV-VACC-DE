using DrWatson
using LibGit2: LibGit2

function apply_patch!(patch, gitpath=projectdir())
    if DrWatson.isdirty(gitpath)
        error("cannot apply patch when HEAD is dirty; please commit or stash changes")
    end

    # Write the patch to file.
    tmpname = tempname() * ".diff"
    open(tmpname, "w") do io
        write(io, patch)
    end

    cmd = `git apply $(tmpname)`
    return cd(gitpath) do
        run(cmd)
    end
end

function checkout!(commitid, gitpath=projectdir())
    if DrWatson.isdirty(gitpath)
        error("cannot apply patch when HEAD is dirty; please commit or stash changes")
    end

    repo = LibGit2.GitRepo(gitpath)
    LibGit2.checkout!(repo, commitid)
end

function checkout!(args::Dict, gitpath=projectdir())
    gitcommit = args["gitcommit"]
    postfix = gitcommit[end - 5:end]
    # DrWatson changed from "_dirty" ↦ "-dirty", so let's stay compatible with both.
    wasdirty = postfix == "_dirty" || postfix == "-dirty"
    # Extract the actual commit.
    gitcommit = !wasdirty ? gitcommit : gitcommit[1:end - 6]

    # Checkout the commit.
    checkout!(gitcommit, gitpath)

    # Apply patch if necessary.
    if wasdirty
        if !haskey(args, "gitpatch")
            @warn "source was dirty but no `gitpatch` in `args`' can't apply patch"
        else
            apply_patch!(args["gitpatch"], gitpath)
        end
    end
end
