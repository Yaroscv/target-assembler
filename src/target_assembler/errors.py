class TargetAssemblerError(RuntimeError):
    """A controlled pipeline failure that should produce a non-zero exit code."""


class ConfigError(TargetAssemblerError):
    """Invalid or inconsistent user configuration."""


class ExternalToolError(TargetAssemblerError):
    """An external bioinformatics program failed."""

