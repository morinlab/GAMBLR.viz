# Pull Request Checklists

**Important:** When opening a pull request, keep only the applicable checklist and delete all other sections.

## Checklist for all PRs

### Required

- [ ] I ensured not to prepend the functions existing in both GAMBLR.data and GAMBLR.results with `<package>::function()` syntax.
- [ ] I ensured not to add GAMBLR.data or GAMBLR.results to the `@import` section of the documentation
- [ ] I tested the new function/functionality in a fresh workspace
- [ ] I added at least one working example that demonstrates the new functionality (if applicable)
- [ ] The test script tools/logExampleOutputs.R ran to completion
- [ ] I included the file GAMBLR_examples_output.log in my pull request and confirm the changes to this file are expected
