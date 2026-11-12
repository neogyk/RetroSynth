"""Language model for the agents: the dspy API model, or the local `claude` CLI on a Claude subscription.

Each agent section in config.ini sets `model_type`:
- `api` (default): `model` is a dspy model string such as anthropic/claude-sonnet-5-5, billed per token
- `subscription`: `model` is ignored and each call runs `claude -p` (prompt on stdin) with the CLI's own login
"""

import json
import subprocess
import tempfile

import dspy
from dspy.lm15 import Message, Response, TextPart, ToolCallPart, ToolResultPart, Usage

API = "api"
SUBSCRIPTION = "subscription"
MODEL_TYPES = (API, SUBSCRIPTION)


def _text(parts):
    """The text of a message's parts. Tool calls and tool results are written out as text, because
    the CLI reads plain text and would otherwise never see its earlier calls or their results."""
    if isinstance(parts, str):
        return parts
    texts, calls = [], []
    for part in parts:
        if isinstance(part, TextPart):
            texts.append(part.text)
        elif isinstance(part, ToolCallPart):
            calls.append({"name": part.name, "args": part.input})
        elif isinstance(part, ToolResultPart):
            label = "error" if part.is_error else "result"
            texts.append(f"<tool_result name={part.name!r} id={part.id!r} {label}>\n{_text(part.content)}\n</tool_result>")
    if calls:
        texts.append(f"<tool_calls>\n{json.dumps({'tool_calls': calls})}\n</tool_calls>")
    return "\n".join(texts)


def _request_text(request):
    """The chat messages as one string of `role: content` lines, as the CLI reads them."""
    return "\n".join(f"{m.role}: {_text(m.parts)}" for m in request.messages)


class ClaudeCLIEngine:
    """dspy engine (DSPy >= 3.5) that answers each request through `claude -p`.

    The dspy system prompt replaces Claude Code's own, the CLI's tools are off (dspy runs the tools),
    and no session is saved. The prompt goes in on stdin: one argument is capped at 128 KB on Linux,
    and the history grows every turn. The CLI runs in the temp directory, so it does not read this
    project's CLAUDE.md or memory."""

    def complete(self, request):
        command = ["claude", "-p", "--tools", "", "--no-session-persistence"]
        if request.system is not None:
            command += ["--system-prompt", _text(request.system)]
        result = subprocess.run(
            command, input=_request_text(request), capture_output=True, text=True, cwd=tempfile.gettempdir(),
        )
        if result.returncode != 0:
            raise RuntimeError(f"claude -p failed (exit {result.returncode}): {result.stderr.strip() or result.stdout.strip()}")
        return Response(
            id=None, model="claude-subscription", finish_reason="stop",
            message=Message.assistant([TextPart(result.stdout.strip())]),
            usage=Usage(input_tokens=0, output_tokens=0, total_tokens=0),
        )


class ClaudeSubscriptionLM(dspy.BaseLM):
    """dspy LM that answers through the `claude` CLI, using the subscription login on this machine."""

    def __init__(self, max_tokens=None, **kwargs):
        super().__init__(model="claude-subscription", max_tokens=max_tokens, engine=ClaudeCLIEngine(), **kwargs)


def make_lm(model, model_type, max_tokens):
    """The dspy LM for one agent, chosen by its config `model_type`."""
    if model_type == SUBSCRIPTION:
        return ClaudeSubscriptionLM(max_tokens=max_tokens)
    if model_type == API:
        return dspy.LM(model, max_tokens=max_tokens)
    raise ValueError(f"unknown model_type {model_type!r}; use one of {list(MODEL_TYPES)}")
