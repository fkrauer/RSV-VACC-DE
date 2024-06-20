# TODO: Convert all this into some function instead.

if !(@isdefined(_args))
    # This allows us to `include` this file while having pre-specified
    # the _unparsed_ arguments `_args`. This means that we can inherit
    # all the defaults, etc. from `parse_args`.
    _args = ARGS
end

args = parse_args(_args, s)
