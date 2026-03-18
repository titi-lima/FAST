### FAST TCP + Vitest Integration

This scaffolding lets you run Vitest tests in a FAST-prioritized order.

### Quick start

1) Initialize at your project root:

```bash
fast-tcp init vitest
```

This creates:
- `.fast/in` and `.fast/out` directories
- `fast-tcp/vitest/` helper scripts
- Adds `npm run test:fast` (if `package.json` exists)

2) Run prioritized tests:

```bash
npm run test:fast
```

What happens:
- Discover test (file, fullName)
- Generate bbox tokens per test case
- Run FAST prioritization via the Python CLI
- Map prioritized IDs back to (file, name)
- Execute Vitest in that order

### Customize

Re-run init with flags to set algorithm and repetitions:

```bash
fast-tcp init vitest --algo FAST-pw --repetitions 3
```

You can also skip package.json modification:

```bash
fast-tcp init vitest --no-scripts
```

The runner lives at `.fast/tools/vitest/run-fast.cjs`. Edit it if you need to pass Vitest flags or batch by file.

### Requirements

- Python 3.9+
- Node 18+
- `fast-tcp` installed in your Python environment

