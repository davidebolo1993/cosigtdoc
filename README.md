# cosigt documentation

Source for the [cosigt documentation](https://davidebolo1993.github.io/cosigtdoc/), built with [VuePress](https://v1.vuepress.vuejs.org/).

## Publishing

Pushing to `master` builds the site and publishes it to the `gh-pages` branch via [`.github/workflows/docs.yml`](.github/workflows/docs.yml).

## Working locally

```bash
# requires node and npm
node --version
npm --version

npm ci

# live preview on http://localhost:8080
npm run docs:dev

# or produce the static site in docs/.vuepress/dist
NODE_OPTIONS=--openssl-legacy-provider npm run docs:build
```

`NODE_OPTIONS=--openssl-legacy-provider` is needed on Node 17 and newer: VuePress 1 builds through webpack 4, whose md4 hashing is rejected by OpenSSL 3.

## Layout

```txt
docs/
├── index.md                            landing page
├── introduction/introduction.md
├── setup/setup.md                      installation, make init/check/run
├── configuration/configuration.md      input tables and every config key
├── workflow/workflow.md                how the pipeline works
├── usecases/usecases.md                worked end-to-end example
├── benchmarking/benchmarking.md        leave-zero-out and leave-all-out
├── contact/contact.md
└── .vuepress/config.js                 nav and sidebar
```

Adding a page means creating the markdown file and registering it in the `sidebar` array in `docs/.vuepress/config.js`.
