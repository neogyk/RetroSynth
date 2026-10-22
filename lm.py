"""Language model for the agents: the dspy API model, or the local `claude` CLI on a Claude subscription.

Each agent section in config.ini sets `model_type`:
- `api` (default): `model` is a dspy model string such as anthropic/claude-sonnet-5-5, billed per token
- `subscription`: `model` is ignored and each call runs `claude -p <prompt>` with the CLI's own login
"""

import subprocess
from types import SimpleNamespace

import dspy

API = "api"
SUBSCRIPTION = "subscription"
MODEL_TYPES = (API, SUBSCRIPTION)


def _messages_text(prompt, messages):
    """The prompt as one string: chat messages become `role: content` lines, as the CLI reads them."""
    if messages:
        return "\n".join(f"{m['role']}: {m['content']}" for m in messages)
    return prompt or ""


def _response(text):
    """An OpenAI-shaped response, which dspy reads for the completion text and usage."""
    message = SimpleNamespace(role="assistant", content=text)
    return SimpleNamespace(
        choices=[SimpleNamespace(index=0, message=message, finish_reason="stop")],
        usage={"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0},
        model="claude-subscription",
    )


class ClaudeSubscriptionLM(dspy.BaseLM):
    """dspy LM that answers through the `claude` CLI, using the subscription login on this machine."""

    def __init__(self, max_tokens=None, **kwargs):
        super().__init__(model="claude-subscription", max_tokens=max_tokens, **kwargs)

    def forward(self, prompt=None, messages=None, **kwargs):
        result = subprocess.run(
            ["claude", "-p", _messages_text(prompt, messages)],
            capture_output=True, text=True, check=True,
        )
        return _response(result.stdout.strip())


def make_lm(model, model_type, max_tokens):
    """The dspy LM for one agent, chosen by its config `model_type`."""
    if model_type == SUBSCRIPTION:
        return ClaudeSubscriptionLM(max_tokens=max_tokens)
    if model_type == API:
        return dspy.LM(model, max_tokens=max_tokens)
    raise ValueError(f"unknown model_type {model_type!r}; use one of {list(MODEL_TYPES)}")
