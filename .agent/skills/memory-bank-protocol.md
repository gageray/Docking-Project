---
name: memory-bank-protocol
description: Protocol for maintaining context across sessions using a memory bank.
---

# Memory Bank Protocol

At the start of every single session, you MUST read the files in the `memory-bank/` folder to boot up your context. At the end of every task, you MUST update `progress.md` and `active-context.md` so you do not forget what you just did.
