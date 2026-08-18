module.exports = {
  title: 'Cosigt documentation',
  base: '/cosigtdoc/',
  theme: 'vt',
  themeConfig: {
    enableDarkMode: true,
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Contact', link: '/contact/contact.md' },
      { text: 'GitHub', link: 'https://github.com/davidebolo1993/cosigt' }
    ],
    sidebar: [
      '/introduction/introduction.md',
      '/setup/setup.md',
      '/configuration/configuration.md',
      '/workflow/workflow.md',
      '/usecases/usecases.md',
      '/benchmarking/benchmarking.md',
      '/contact/contact.md'
    ],
  }
}
