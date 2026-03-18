## FAST TCP Ant Plugin Usage

This example ships a ready-to-use Ant integration for FAST TCP that makes it simple to prioritize and run your JUnit 5 tests in an optimized order.

### Prerequisites

- Java 21 (JDK 21)
- Apache Ant 1.10+
- Python 3.8+
- FAST TCP installed: `pip install fast-tcp`

### Initialize once (opinionated, macros-only)

From the repo root:

```bash
cd examples/java/ant
fast-tcp init ant
```

This will:

- Ensure JUnit Console is available under `lib/`
- Compile the project
- Generate black-box inputs from `src/test/java` into `.fast/in/<project>-bbox.txt`
- Run FAST TCP prioritization (`FAST-pw` by default, 3 repetitions)
- Map prioritized IDs back to JUnit method selectors
- Create `.fast/fast-tcp.xml` pinned to your current Python interpreter
- Inject `<import file=".fast/fast-tcp.xml"/>` and a `fast-tcp` target into `build.xml`
- Execute JUnit Console in that prioritized order

### Recommended: import macros in your build.xml

Add this to the top of your `build.xml` so you can call a single target:

```xml
<property name="fast.tcp.macros" value="${user.home}/.fast-tcp/fast-tcp.xml"/>
<target name="fast-tcp-download-macros">
    <mkdir dir="${user.home}/.fast-tcp"/>
    <get src="https://raw.githubusercontent.com/icse18-FAST/FAST/master/fast_tcp/integrations/ant/fast-tcp.xml" dest="${user.home}/.fast-tcp/fast-tcp.xml"/>
    <echo>FAST TCP macros downloaded to ${fast.tcp.macros}</echo>
</target>

<import file="${fast.tcp.macros}"/>

<target name="fast-tcp" depends="fast-tcp-download-macros">
    <fast-tcp-prioritize-and-run projectDir="${basedir}" algo="FAST-pw" repetitions="3"/>
</target>
```

Now simply run:

```bash
ant fast-tcp
```

### Configurable properties

- During init: `--project-dir`, `--algo`, and `--repetitions`
- After init: adjust the generated `fast-tcp` target in `build.xml` if desired

### Project layout assumptions

- Tests under `src/test/java`
- JUnit 5 Console used to run tests (downloaded automatically to `lib/`)

### Outputs

- `.fast/in/<project>-bbox.txt` — tokenized black-box input
- `.fast/in/selectors.txt` — discovered JUnit5 method selectors
- `.fast/in/prioritized-selectors.txt` — prioritized JUnit selectors
- `.fast/out/<dataset>/prioritized/*` — FAST TCP result files (`.pickle`, `.tsv`)

