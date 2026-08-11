from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, TextIO

from .errors import ExternalToolError


class Runner:
    def __init__(self, logger: logging.Logger, dry_run: bool = False):
        self.logger = logger
        self.dry_run = dry_run

    def require(self, programs: Iterable[str]) -> None:
        missing = sorted({program for program in programs if shutil.which(program) is None})
        if missing:
            raise ExternalToolError("Missing required programs: " + ", ".join(missing))

    def run(
        self,
        command: list[str],
        *,
        cwd: str | Path | None = None,
        stdout_path: str | Path | None = None,
        stderr_path: str | Path | None = None,
        check: bool = True,
    ) -> subprocess.CompletedProcess:
        display = shlex.join(str(part) for part in command)
        self.logger.info("RUN %s", display)
        if self.dry_run:
            return subprocess.CompletedProcess(command, 0)
        stdout_handle: TextIO | None = None
        stderr_handle: TextIO | None = None
        try:
            if stdout_path:
                stdout_handle = open(stdout_path, "w", encoding="utf-8")
            if stderr_path:
                stderr_handle = open(stderr_path, "w", encoding="utf-8")
            result = subprocess.run(
                command,
                cwd=cwd,
                stdout=stdout_handle,
                stderr=stderr_handle,
                text=True,
                check=False,
            )
        finally:
            if stdout_handle:
                stdout_handle.close()
            if stderr_handle:
                stderr_handle.close()
        if check and result.returncode != 0:
            tail = ""
            if stderr_path and Path(stderr_path).exists():
                lines = Path(stderr_path).read_text(encoding="utf-8", errors="replace").splitlines()
                tail = "\n" + "\n".join(lines[-25:])
            raise ExternalToolError(f"Command failed ({result.returncode}): {display}{tail}")
        return result

    def popen_stdout(
        self,
        command: list[str],
        *,
        cwd: str | Path | None = None,
        stderr_path: str | Path,
    ) -> tuple[subprocess.Popen, TextIO]:
        display = shlex.join(str(part) for part in command)
        self.logger.info("RUN/STREAM %s", display)
        if self.dry_run:
            raise ExternalToolError("Streaming commands are unavailable during dry-run")
        stderr_handle = open(stderr_path, "w", encoding="utf-8")
        process = subprocess.Popen(
            command,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        return process, stderr_handle


def check_stream(process: subprocess.Popen, stderr_handle: TextIO, description: str) -> None:
    return_code = process.wait()
    stderr_path = Path(stderr_handle.name)
    stderr_handle.close()
    if return_code != 0:
        lines = stderr_path.read_text(encoding="utf-8", errors="replace").splitlines()
        raise ExternalToolError(
            f"{description} failed ({return_code})\n" + "\n".join(lines[-25:])
        )

