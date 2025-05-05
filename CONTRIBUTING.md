# Contributing to Orbweaver

Thank you for your interest in contributing to Orbweaver! This document provides guidelines and instructions for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Environment](#development-environment)
- [Development Workflow](#development-workflow)
- [Code Style Guidelines](#code-style-guidelines)
- [Testing Guidelines](#testing-guidelines)
- [Documentation Guidelines](#documentation-guidelines)
- [Pull Request Process](#pull-request-process)
- [Issue Reporting](#issue-reporting)
- [Feature Requests](#feature-requests)
- [Project Structure](#project-structure)
- [Communication](#communication)

## Code of Conduct

By participating in this project, you agree to maintain a respectful and inclusive environment for all contributors. Please be considerate in your communications and respect different perspectives and experiences.

## Getting Started

1. **Fork the Repository**
   - Click the "Fork" button at the top right of the repository page.

2. **Clone Your Fork**
   ```bash
   git clone https://github.com/YOUR-USERNAME/orbweaver.git
   cd orbweaver
   ```

3. **Add Upstream Remote**
   ```bash
   git remote add upstream https://github.com/ORIGINAL-OWNER/orbweaver.git
   ```

4. **Build the Project**
   ```bash
   cargo build
   ```

5. **Run Tests**
   ```bash
   cargo test
   ```

## Development Environment

### Recommended Tools

- **Rust**: Latest stable version (1.70.0+)
- **IDE**: VS Code with rust-analyzer
- **Linting**: `cargo clippy`
- **Formatting**: `cargo fmt`
- **Documentation**: `cargo doc --open`
- **Optional**: Graphviz for visualization

### Development Dependencies

- Install Rust using [rustup](https://rustup.rs/)
- Install development tools:
  ```bash
  rustup component add clippy rustfmt
  cargo install cargo-watch cargo-audit
  ```

## Development Workflow

1. **Create a Branch**
   - Create a branch for your work:
     ```bash
     git checkout -b feature/your-feature-name
     ```
   - Use prefixes like `feature/`, `bugfix/`, `docs/`, etc.

2. **Make Changes**
   - Write your code following our code style guidelines
   - Keep your changes focused and related to the issue/feature

3. **Test Your Changes**
   - Add or update tests
   - Run existing tests:
     ```bash
     cargo test
     ```
   - Run linting:
     ```bash
     cargo clippy -- -D warnings
     ```
   - Format your code:
     ```bash
     cargo fmt
     ```

4. **Commit Your Changes**
   - Use descriptive commit messages
   - Reference issue numbers in commit messages
   - Example:
     ```
     Add reverse complement support for minimizers
     
     - Implement strand-aware minimizer generation
     - Add tests for reverse complement cases
     - Update documentation
     
     Fixes #42
     ```

5. **Keep Updated with Upstream**
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

6. **Push to Your Fork**
   ```bash
   git push origin feature/your-feature-name
   ```

7. **Create a Pull Request**
   - See [Pull Request Process](#pull-request-process) below

8. **Address Review Feedback**
   - Make requested changes
   - Commit and push updates
   - Discuss any disagreements respectfully

## Code Style Guidelines

Orbweaver follows Rust standard formatting and style conventions:

1. **Formatting**:
   - Use `cargo fmt` to automatically format your code
   - Do not disable formatting for blocks unless absolutely necessary

2. **Naming Conventions**:
   - Use `snake_case` for variables, functions, and modules
   - Use `CamelCase` for types and traits
   - Use `SCREAMING_SNAKE_CASE` for constants
   - Prefix unsafe functions with `unsafe_`

3. **Comments**:
   - Use `///` for documentation comments on public items
   - Document all public interfaces
   - Include examples in documentation for complex functions
   - Explain *why* rather than *what* for internal code comments

4. **Error Handling**:
   - Use `anyhow::Result` for function returns that can fail
   - Add context to errors using `.context()` or `.with_context()`
   - Avoid panics in library code

5. **Clippy**:
   - Follow clippy suggestions
   - Fix clippy warnings before submitting PRs

## Testing Guidelines

1. **Test Coverage**:
   - Add tests for all new functionality
   - Aim for at least 80% coverage for new code
   - Include both positive and negative test cases

2. **Test Types**:
   - **Unit Tests**: Test individual functions and methods
   - **Integration Tests**: Test workflow across modules
   - **Documentation Tests**: Ensure examples in docs work

3. **Test Organization**:
   - Unit tests should be in the same file as the code being tested, inside a `#[cfg(test)]` module
   - Integration tests go in the `tests/` directory
   - Use descriptive test names that explain the scenario and expected outcome

4. **Running Tests**:
   ```bash
   # Run all tests
   cargo test
   
   # Run specific test
   cargo test test_name
   
   # Run tests with output
   cargo test -- --nocapture
   ```

## Documentation Guidelines

1. **Code Documentation**:
   - Document all public items with `///` comments
   - Include parameters, return values, and examples
   - Explain the purpose and usage of the item
   - Mention any edge cases or performance considerations

2. **Module Documentation**:
   - Add `//!` comments at the top of each module file
   - Explain the purpose and components of the module
   - Show a simple usage example if applicable

3. **Project Documentation**:
   - Update README.md with significant changes
   - Update architecture.md if you change the design
   - Add to output_formats.md if you add or modify outputs
   - Create new docs for major features

## Pull Request Process

1. **Create a Pull Request**:
   - Go to the original repository
   - Click "New Pull Request"
   - Select "compare across forks"
   - Select your fork and branch

2. **PR Description**:
   - Use the PR template if provided
   - Describe what changes you've made
   - Explain why they are needed
   - Reference any related issues
   - Include any testing instructions
   - Link to any relevant documentation

3. **Review Process**:
   - Maintainers will review your code
   - Address any feedback
   - Make requested changes
   - Respond to comments

4. **CI Checks**:
   - All automated tests must pass
   - No clippy warnings
   - Code must be correctly formatted

5. **Merge**:
   - A maintainer will merge your PR once approved
   - Typically, PRs will be squash-merged

## Issue Reporting

1. **Search First**:
   - Check if the issue has already been reported

2. **Issue Template**:
   - Use the appropriate issue template if available
   - If not, include:
     - Description of the problem
     - Steps to reproduce
     - Expected behavior
     - Actual behavior
     - Environment details (OS, Rust version, etc.)
     - Screenshots or logs if applicable

3. **Labels**:
   - Use appropriate labels if you have permission
   - Common labels: bug, enhancement, documentation, help wanted

## Feature Requests

1. **Describe the Need**:
   - Explain what problem the feature would solve
   - Describe why it's important

2. **Outline the Solution**:
   - Provide a high-level description of how it might work
   - Note any alternatives you've considered

3. **Implementation Ideas** (optional):
   - Any thoughts on how it could be implemented
   - Potential challenges or trade-offs

## Project Structure

To understand the project structure better, refer to:
- [docs/architecture.md](docs/architecture.md) - Overall architecture
- [docs/cli_usage.md](docs/cli_usage.md) - CLI usage information
- [docs/output_formats.md](docs/output_formats.md) - Output formats

Key directories:
- `src/` - Source code
- `docs/` - Documentation
- `tests/` - Integration tests
- `scripts/` - Utility scripts

## Communication

- **Issues**: Use for bug reports and feature requests
- **Pull Requests**: Use for code contributions
- **Discussions**: Use for general questions and discussions

---

Thank you for contributing to Orbweaver! Your efforts help improve genomic sequence analysis for everyone. 