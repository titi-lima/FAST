# Publishing to PyPI

This document explains how to set up automatic publishing to PyPI using GitHub Actions.

## Prerequisites

### 1. PyPI Account Setup

1. **Create PyPI Account**: Register at [pypi.org](https://pypi.org/account/register/)
2. **Create Test PyPI Account**: Register at [test.pypi.org](https://test.pypi.org/account/register/)
3. **Enable 2FA**: Enable two-factor authentication on both accounts

### 2. Generate API Tokens

#### For PyPI (Production)
1. Go to [PyPI Account Settings](https://pypi.org/manage/account/)
2. Scroll to "API tokens" section
3. Click "Add API token"
4. Set token name: `fast-tcp-github-actions`
5. Set scope: `Entire account` (or specific to `fast-tcp` project once created)
6. **Copy the token** (starts with `pypi-`)

#### For Test PyPI
1. Go to [Test PyPI Account Settings](https://test.pypi.org/manage/account/)
2. Follow the same steps as above
3. Set token name: `fast-tcp-github-actions-test`
4. **Copy the token** (starts with `pypi-`)

### 3. Configure GitHub Secrets

1. Go to your GitHub repository
2. Navigate to **Settings** → **Secrets and variables** → **Actions**
3. Add the following secrets:

| Secret Name | Value | Description |
|-------------|--------|-------------|
| `PYPI_API_TOKEN` | `pypi-AgEIcHl...` | Production PyPI API token |
| `TEST_PYPI_API_TOKEN` | `pypi-AgEIcHl...` | Test PyPI API token |

### 4. Set Up GitHub Environments (Optional but Recommended)

1. Go to **Settings** → **Environments**
2. Create environment: `pypi`
   - Add protection rule: Required reviewers (yourself)
   - Add environment secret: `PYPI_API_TOKEN`
3. Create environment: `testpypi`
   - Add environment secret: `TEST_PYPI_API_TOKEN`

## Publishing Workflow

The GitHub Actions workflow (`publish.yml`) supports multiple publishing scenarios:

### 1. Automatic Publishing on Release

**For Production PyPI:**
```bash
# Create a new release
git tag v1.0.0
git push origin v1.0.0

# Or create release via GitHub UI
# - Go to Releases → Create a new release
# - Tag: v1.0.0
# - Title: v1.0.0
# - Description: Release notes
# - Publish release
```

**For Test PyPI (Pre-releases):**
```bash
# Create a pre-release
git tag v1.0.0-rc1
git push origin v1.0.0-rc1

# In GitHub UI, mark as "pre-release"
```

### 2. Manual Publishing

Trigger the workflow manually:
1. Go to **Actions** tab
2. Select "Publish to PyPI" workflow
3. Click "Run workflow"
4. This will publish to Test PyPI only

## Version Management

### Update Version Numbers

Before creating a release, update version numbers in:

1. **fast_tcp/__init__.py**:
```python
__version__ = "1.0.1"
```

2. **pyproject.toml**:
```toml
version = "1.0.1"
```

3. **setup.py**:
```python
version="1.0.1",
```

### Version Scheme

Follow [Semantic Versioning](https://semver.org/):
- `MAJOR.MINOR.PATCH` (e.g., `1.0.0`)
- `MAJOR.MINOR.PATCH-rc.N` for release candidates (e.g., `1.0.0-rc.1`)
- `MAJOR.MINOR.PATCH-alpha.N` for alpha versions
- `MAJOR.MINOR.PATCH-beta.N` for beta versions

## Workflow Details

The publish workflow includes three jobs:

### 1. Build (`build`)
- Builds source distribution and wheel
- Runs on Ubuntu with Python 3.x
- Uploads build artifacts

### 2. Test Build (`test-build`)
- Tests installation across Python 3.8-3.12
- Verifies imports work correctly
- Tests CLI functionality

### 3. Publish
- **Test PyPI**: For pre-releases and manual triggers
- **Production PyPI**: For stable releases only

## Troubleshooting

### Common Issues

**1. "File already exists" error:**
- PyPI doesn't allow re-uploading the same version
- Increment version number and try again

**2. "Invalid credentials" error:**
- Check API token is correct
- Ensure token has proper scope
- Verify secret name matches workflow

**3. "Package name already taken":**
- Choose a different package name
- Update `name` in `setup.py` and `pyproject.toml`

**4. Build failures:**
- Check dependencies in `requirements.txt`
- Ensure all necessary files are included in `MANIFEST.in`

### Testing Before Release

**Test locally:**
```bash
# Build packages
python -m build

# Test installation
pip install dist/fast_tcp-*.whl

# Test functionality
fast-tcp --help
python -c "import fast_tcp; print(fast_tcp.__version__)"
```

**Test on Test PyPI:**
```bash
# Upload to Test PyPI manually
python -m twine upload --repository testpypi dist/*

# Install from Test PyPI
pip install --index-url https://test.pypi.org/simple/ fast-tcp
```

## Security Best Practices

1. **Use API tokens**, not username/password
2. **Scope tokens** to specific projects when possible
3. **Use GitHub environments** for additional protection
4. **Enable 2FA** on PyPI accounts
5. **Regular token rotation** (every 6-12 months)
6. **Monitor releases** for unauthorized uploads

## Release Checklist

Before creating a release:

- [ ] Update version numbers in all files
- [ ] Update CHANGELOG.md (if you have one)
- [ ] Test package locally
- [ ] Run tests
- [ ] Create and test pre-release if needed
- [ ] Create GitHub release
- [ ] Verify PyPI upload was successful
- [ ] Test installation from PyPI

## Example Release Process

```bash
# 1. Update version
sed -i 's/__version__ = ".*"/__version__ = "1.0.1"/' fast_tcp/__init__.py
sed -i 's/version = ".*"/version = "1.0.1"/' pyproject.toml
sed -i 's/version=".*"/version="1.0.1"/' setup.py

# 2. Commit changes
git add .
git commit -m "Bump version to 1.0.1"
git push

# 3. Create and push tag
git tag v1.0.1
git push origin v1.0.1

# 4. Create GitHub release (via UI or CLI)
gh release create v1.0.1 --title "v1.0.1" --notes "Release notes here"
```

The GitHub Actions workflow will automatically:
1. Build the package
2. Test across multiple Python versions
3. Publish to PyPI
4. Create release artifacts

Your package will be available at: https://pypi.org/project/fast-tcp/ 