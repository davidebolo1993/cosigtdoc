module.exports = {
  title: 'Template',
  base: '/vuepress_doc_template/',
  theme: 'vt',
  themeConfig: {
    enableDarkMode: true,
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Contact', link: '/contact/contact.md' },
      { text: 'GitHub', link: 'https://github.com/davidebolo1993/vuepress_doc_template' }
    ],
    sidebar: [
      '/introduction/introduction.md',
      '/setup/setup.md',
      '/contact/contact.md'
    ],
  }
}

