var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _toConsumableArray(arr) { if (Array.isArray(arr)) { for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) { arr2[i] = arr[i]; } return arr2; } else { return Array.from(arr); } }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

function Header(props) {
  return React.createElement(
    "nav",
    { className: "navbar navbar-expand-md sticky-top navbar-light", style: { backgroundColor: "#e3f2fd" } },
    React.createElement(
      "div",
      { className: "navbar-collapse collapse w-100 order-1 order-md-0 dual-collapse2" },
      React.createElement(
        "ul",
        { className: "navbar-nav mr-auto" },
        React.createElement("img", { src: "img/bbi_icon.png", height: "70", className: "d-inline-block align-top", alt: "" })
      )
    ),
    React.createElement(
      "div",
      { className: "mx-auto order-0" },
      React.createElement(
        "a",
        { className: "navbar-brand mx-auto", href: "#" },
        "Experiment ",
        props.run_name,
        " QC Dashboard"
      )
    ),
    React.createElement("div", { className: "navbar-collapse collapse w-100 order-3 dual-collapse2" })
  );
}

/*
** Sample-specific display on right side of page.
*/
function Sample(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  var barnyardRegex = new RegExp("^[bB]arnyard");
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: safe_name, role: "tabpanel", "aria-labelledby": safe_name },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
      React.createElement(
        "h1",
        { className: "h3", id: "lig-name" },
        props.sample_id
      )
    ),
    React.createElement(
      "nav",
      null,
      React.createElement(
        "div",
        { className: "nav nav-tabs", id: "nav" + safe_name + "-tab", role: "tablist" },
        barnyardRegex.test(props.sample_id) && React.createElement(
          "a",
          {
            className: "nav-item nav-link active",
            id: "nav" + safe_name + "-barn-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-barn",
            role: "tab", "aria-controls": "nav" + safe_name + "-barn",
            "aria-selected": "true" },
          "Barnyard"
        ),
        barnyardRegex.test(props.sample_id) ? React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-read_distribution-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-read_distribution",
            role: "tab", "aria-controls": "nav" + safe_name + "-read_distribution",
            "aria-selected": "true" },
          "Read distributions"
        ) : React.createElement(
          "a",
          {
            className: "nav-item nav-link active",
            id: "nav" + safe_name + "-read_distribution-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-read_distribution",
            role: "tab", "aria-controls": "nav" + safe_name + "-read_distribution",
            "aria-selected": "true" },
          "Read distributions"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-scrub-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-scrub",
            role: "tab", "aria-controls": "nav" + safe_name + "-scrub",
            "aria-selected": "false" },
          "Scrublet"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-frit-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-frit",
            role: "tab", "aria-controls": "nav" + safe_name + "-frit",
            "aria-selected": "false" },
          "FRIT"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-umap-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-umap",
            role: "tab", "aria-controls": "nav" + safe_name + "-umap",
            "aria-selected": "false" },
          "UMAP"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-fragmentdist-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-fragmentdist",
            role: "tab", "aria-controls": "nav" + safe_name + "-fragmentdist",
            "aria-selected": "false" },
          "Estimated fragment size"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-frip-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-frip",
            role: "tab", "aria-controls": "nav" + safe_name + "-frip",
            "aria-selected": "false" },
          "FRIP"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-insertsize-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-insertsize",
            role: "tab", "aria-controls": "nav" + safe_name + "-insertsize",
            "aria-selected": "false" },
          "Insert size"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-tssenrichment-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-tssenrichment",
            role: "tab", "aria-controls": "nav" + safe_name + "-tssenrichment",
            "aria-selected": "false" },
          "TSS Enrichment"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-stats-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-stats",
            role: "tab", "aria-controls": "nav" + safe_name + "-stats",
            "aria-selected": "false" },
          "Sample Stats"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-readmetrics-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-readmetrics",
            role: "tab", "aria-controls": "nav" + safe_name + "-readmetrics",
            "aria-selected": "false" },
          "Read Metrics"
        ),
        React.createElement(
          "a",
          {
            className: "nav-item nav-link",
            id: "nav" + safe_name + "-fulllog-tab",
            "data-toggle": "tab", href: "#nav" + safe_name + "-fulllog",
            role: "tab", "aria-controls": "nav" + safe_name + "-fulllog",
            "aria-selected": "false" },
          "Full Log"
        )
      )
    ),
    React.createElement(
      "div",
      { className: "tab-content", id: "nav-tabContent" },
      barnyardRegex.test(props.sample_id) && React.createElement(BarnyardPane, { sample_id: props.sample_id, className: "tab-pane fade show active" }),
      barnyardRegex.test(props.sample_id) ? React.createElement(ReadDistributionsPane, { sample_id: props.sample_id, className: "tab-pane fade" }) : React.createElement(ReadDistributionsPane, { sample_id: props.sample_id, className: "tab-pane fade show active" }),
      React.createElement(ScrubPane, { sample_id: props.sample_id, sample_stats: run_data.sample_stats }),
      React.createElement(FritPane, { sample_id: props.sample_id }),
      React.createElement(UmapPane, { sample_id: props.sample_id }),
      React.createElement(FragmentDistPane, { sample_id: props.sample_id }),
      React.createElement(FripPane, { sample_id: props.sample_id }),
      React.createElement(InsertSizePane, { sample_id: props.sample_id }),
      React.createElement(TssEnrichmentPane, { sample_id: props.sample_id }),
      React.createElement(StatsPane, { sample_id: props.sample_id, sample_stats: run_data.sample_stats }),
      React.createElement(ReadMetricsPane, { sample_id: props.sample_id, log: log_data[props.sample_id] }),
      React.createElement(FullLogPane, { sample_id: props.sample_id, log: full_log_data[props.sample_id] })
    )
  );
}

/*
** Display a image in a window pane.
*/
function Pane(props) {
  return React.createElement(
    "div",
    { className: props.className, id: props.id, role: "tabpanel", "aria-labelledby": props.tag },
    React.createElement(
      "p",
      null,
      props.text.map(function (text, index) {
        return typeof text == "string" ? React.createElement(
          "span",
          { key: index },
          text
        ) : React.createElement(
          "a",
          { key: index, href: text.link },
          text.label
        );
      })
    ),
    React.createElement("img", { src: props.plot, className: "rounded mx-auto d-block", alt: "...", style: { maxHeight: "50vh", width: "auto" } })
  );
}

/*
** Display multiple images in a window pane.
*/
function PaneMultiImage(props) {
  return React.createElement(
    "div",
    { className: props.className, id: props.id, role: "tabpanel", "aria-labelledby": props.tag },
    React.createElement(
      "p",
      null,
      props.text.map(function (text, index) {
        return typeof text == "string" ? React.createElement(
          "span",
          { key: index },
          text
        ) : React.createElement(
          "a",
          { key: index, href: text.link },
          text.label
        );
      })
    ),
    props.plot_list.map(function (item) {
      return React.createElement("img", { src: item, key: item, className: "rounded mx-auto d-block", alt: "...", style: { height: "auto", width: "auto" } });
    })
  );
}

function TitleRow(props) {
  return React.createElement(
    "th",
    { scope: "col" },
    props.samp
  );
}

function RegRow(props) {
  return React.createElement(
    "td",
    null,
    props.val
  );
}

/*
** Display sample-specific summary statistics table in a sample (tab) pane.
*/
function StatsPane(props) {
  var sample_stat = props.sample_stats[props.sample_id];
  var stats_list = ["Total Reads", "Total merged peaks", "Median reads per cell", "Median per cell FRIP", "Median per cell FRIT", "Median Duplication rate"];
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "nav" + safe_name + "-stats", role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-stats-tab" },
    React.createElement(
      "table",
      { className: "table table-hover" },
      React.createElement(
        "thead",
        null,
        React.createElement(
          "tr",
          null,
          React.createElement("th", { scope: "col" }),
          stats_list.map(function (item, index) {
            return React.createElement(TitleRow, { key: index, samp: item });
          })
        )
      ),
      React.createElement(
        "tbody",
        null,
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row" },
            props.sample_id
          ),
          React.createElement(RegRow, { val: sample_stat.Total_reads }),
          React.createElement(RegRow, { val: sample_stat.Total_merged_peaks }),
          React.createElement(RegRow, { val: sample_stat.Median_reads_per_cell }),
          React.createElement(RegRow, { val: sample_stat.Median_per_cell_frip }),
          React.createElement(RegRow, { val: sample_stat.Median_per_cell_frit }),
          React.createElement(RegRow, { val: sample_stat.Median_duplication_rate })
        )
      )
    )
  );
}

function CodeChunk(props) {
  return React.createElement(
    "pre",
    { style: { paddingLeft: '20px', paddingRight: '20px' } },
    React.createElement(
      "code",
      null,
      '\n' + props.text + '\n\n'
    )
  );
}

function ReadMetricsPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "nav" + safe_name + "-readmetrics", role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-readmetrics-tab" },
    React.createElement(CodeChunk, { text: props.log })
  );
}

function FullLogPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "nav" + safe_name + "-fulllog", role: "tabpanel", "aria-labelledby": "nav" + safe_name + "-fulllog-tab" },
    React.createElement(CodeChunk, { text: props.log })
  );
}

function BarnyardPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(Pane, {
    className: props.className,
    id: "nav" + safe_name + "-barn",
    tag: "nav" + safe_name + "-barn-tab",
    text: [''],
    plot: "img/" + props.sample_id + ".png"
  });
}
function FritPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(Pane, {
    className: "tab-pane fade",
    id: "nav" + safe_name + "-frit",
    tag: "nav" + safe_name + "-frit-tab",
    text: [''],
    plot: "img/" + props.sample_id + "-frit.png"
  });
}
function ScrubPane(props) {
  var sample_stat = props.sample_stats[props.sample_id];
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(Pane, {
    className: "tab-pane fade",
    id: "nav" + safe_name + "-scrub",
    tag: "nav" + safe_name + "-scrub-tab",
    text: ['Doublet count: ' + sample_stat.Doublet_Number + "\n\nDoublet rate: " + sample_stat.Doublet_Percent],
    plot: "img/" + props.sample_id + "-scrublet_hist.png"
  });
}
function ReadDistributionsPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(PaneMultiImage, {
    className: props.className,
    id: "nav" + safe_name + "-read_distribution",
    tag: "nav" + safe_name + "-read_distribution-tab",
    text: [''],
    plot_list: ["img/" + props.sample_id + "-knee_plot.png", "img/" + props.sample_id + "-umi_per_cell.png"]
  });
}
function UmapPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(Pane, {
    className: "tab-pane fade",
    id: "nav" + safe_name + "-umap",
    tag: "nav" + safe_name + "-umap-tab",
    text: [''],
    plot: "img/" + props.sample_id + "-umap_plot.png"
  });
}
function FragmentDistPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(Pane, {
    className: "tab-pane fade",
    id: "nav" + safe_name + "-fragmentdist",
    tag: "nav" + safe_name + "-fragmentdist-tab",
    text: [''],
    plot: "img/" + props.sample_id + "-estimated_frag_dist.png"
  });
}
function FripPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(Pane, {
    className: "tab-pane fade",
    id: "nav" + safe_name + "-frip",
    tag: "nav" + safe_name + "-frip-tab",
    text: [''],
    plot: "img/" + props.sample_id + "-frip.png"
  });
}
function InsertSizePane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(Pane, {
    className: "tab-pane fade",
    id: "nav" + safe_name + "-insertsize",
    tag: "nav" + safe_name + "-insertsize-tab",
    text: [''],
    plot: "img/" + props.sample_id + "-insert_size_dist.png"
  });
}
function TssEnrichmentPane(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(Pane, {
    className: "tab-pane fade",
    id: "nav" + safe_name + "-tssenrichment",
    tag: "nav" + safe_name + "-tssenrichment-tab",
    text: [''],
    plot: "img/" + props.sample_id + "-tss_enrichment.png"
  });
}

function SamplePill(props) {
  var safe_name = "hp" + props.sample_id.replace(/[.]/g, "-");
  return React.createElement(
    "a",
    { className: "nav-link", id: safe_name + "-tab", "data-toggle": "pill", href: "#" + safe_name, role: "tab",
      "aria-controls": safe_name, "aria-selected": "false" },
    props.sample_id
  );
}

//https://www.florin-pop.com/blog/2019/07/sort-table-data-with-react/

var tableData = Object.values(run_data['sample_stats']);

var sortTypes = {
  total_reads_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return a.Total_reads - b.Total_reads;
    }
  },
  total_reads_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return b.Total_reads - a.Total_reads;
    }
  },
  total_peaks_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return a.Total_merged_peaks - b.Total_merged_peaks;
    }
  },
  total_peaks_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return b.Total_merged_peaks - a.Total_merged_peaks;
    }
  },
  sample_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return ('' + a.Sample).localeCompare(b.Sample);
    }
  },
  sample_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return ('' + b.Sample).localeCompare(a.Sample);
    }
  },
  median_dup_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return parseFloat(a.Median_duplication_rate) - parseFloat(b.Median_duplication_rate);
    }
  },
  median_dup_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return parseFloat(b.Median_duplication_rate) - parseFloat(a.Median_duplication_rate);
    }
  },
  median_reads_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return a.Median_reads_per_cell - b.Median_reads_per_cell;
    }
  },
  median_reads_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return b.Median_reads_per_cell - a.Median_reads_per_cell;
    }
  },
  median_frip_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return parseFloat(b.Median_per_cell_frip) - parseFloat(a.Median_per_cell_frip);
    }
  },
  median_frip_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return parseFloat(a.Median_per_cell_frip) - parseFloat(b.Median_per_cell_frip);
    }
  },
  median_frit_up: {
    class: 'sort-up',
    fn: function fn(a, b) {
      return parseFloat(b.Median_per_cell_frit) - parseFloat(a.Median_per_cell_frit);
    }
  },
  median_frit_down: {
    class: 'sort-down',
    fn: function fn(a, b) {
      return parseFloat(a.Median_per_cell_frit) - parseFloat(b.Median_per_cell_frit);
    }
  },
  default: {
    class: 'sort',
    fn: function fn(a, b) {
      return a;
    }
  }
};

/*
** Looks like mechanisms for sorting summay table by user-selected column contents.
*/

var Table = function (_React$Component) {
  _inherits(Table, _React$Component);

  function Table() {
    var _ref;

    var _temp, _this, _ret;

    _classCallCheck(this, Table);

    for (var _len = arguments.length, args = Array(_len), _key = 0; _key < _len; _key++) {
      args[_key] = arguments[_key];
    }

    return _ret = (_temp = (_this = _possibleConstructorReturn(this, (_ref = Table.__proto__ || Object.getPrototypeOf(Table)).call.apply(_ref, [this].concat(args))), _this), _this.state = {
      currentSort: 'sample_up'
    }, _this.onSortTotalReads = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'total_reads_down') nextSort = 'total_reads_up';else if (currentSort === 'total_reads_up') nextSort = 'total_reads_down';else nextSort = 'total_reads_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortSample = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'sample_down') nextSort = 'sample_up';else if (currentSort === 'sample_up') nextSort = 'sample_down';else nextSort = 'sample_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortTotalPeaks = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'total_peaks_down') nextSort = 'total_peaks_up';else if (currentSort === 'total_peaks_up') nextSort = 'total_peaks_down';else nextSort = 'total_peaks_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortMedianDuplication = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'median_dup_down') nextSort = 'median_dup_up';else if (currentSort === 'median_dup_up') nextSort = 'median_dup_down';else nextSort = 'median_dup_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortMedianReads = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'median_reads_down') nextSort = 'median_reads_up';else if (currentSort === 'median_reads_up') nextSort = 'median_reads_down';else nextSort = 'median_reads_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortMedianFrip = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'median_frip_down') nextSort = 'median_frip_up';else if (currentSort === 'median_frip_up') nextSort = 'median_frip_down';else nextSort = 'median_frip_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _this.onSortMedianFrit = function () {
      var currentSort = _this.state.currentSort;

      var nextSort = void 0;

      if (currentSort === 'median_frit_down') nextSort = 'median_frit_up';else if (currentSort === 'median_frit_up') nextSort = 'median_frit_down';else nextSort = 'median_frit_up';

      _this.setState({
        currentSort: nextSort
      });
    }, _temp), _possibleConstructorReturn(_this, _ret);
  }

  // declaring the default state


  // method called every time the sort button is clicked
  // it will change the currentSort value to the next one


  _createClass(Table, [{
    key: "render",


    /*
    ** Show the overall summary table (selected on the left side of page).
    */
    value: function render() {
      var data = this.props.data;
      var currentSort = this.state.currentSort;

      return React.createElement(
        "div",
        { className: "tab-pane fade show active", id: "summary", role: "tabpanel" },
        React.createElement(
          "div",
          { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
          React.createElement(
            "h1",
            { className: "h3", id: "lig-name" },
            "Summary Table"
          )
        ),
        data.length > 0 && React.createElement(
          "table",
          { className: "table table-hover table-responsive summary-table" },
          React.createElement(
            "thead",
            null,
            React.createElement(
              "tr",
              null,
              React.createElement(
                "th",
                null,
                "Sample",
                React.createElement(
                  "button",
                  { onClick: this.onSortSample, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Total reads",
                React.createElement(
                  "button",
                  { onClick: this.onSortTotalReads, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Total merged peaks",
                React.createElement(
                  "button",
                  { onClick: this.onSortTotalPeaks, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Median reads per cell",
                React.createElement(
                  "button",
                  { onClick: this.onSortMedianReads, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Median per cell FRIP",
                React.createElement(
                  "button",
                  { onClick: this.onSortMedianFrip, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Median per cell FRIT",
                React.createElement(
                  "button",
                  { onClick: this.onSortMedianFrit, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              ),
              React.createElement(
                "th",
                null,
                "Median duplication rate",
                React.createElement(
                  "button",
                  { onClick: this.onSortMedianDuplication, className: "sort_button" },
                  React.createElement("i", { className: "fas fa-sort" })
                )
              )
            )
          ),
          React.createElement(
            "tbody",
            null,
            [].concat(_toConsumableArray(data)).sort(sortTypes[currentSort].fn).map(function (p, index) {
              return React.createElement(
                "tr",
                { key: index },
                React.createElement(
                  "td",
                  null,
                  p.Sample
                ),
                React.createElement(
                  "td",
                  null,
                  p.Total_reads
                ),
                React.createElement(
                  "td",
                  null,
                  p.Total_merged_peaks
                ),
                React.createElement(
                  "td",
                  null,
                  p.Median_reads_per_cell
                ),
                React.createElement(
                  "td",
                  null,
                  p.Median_per_cell_frip
                ),
                React.createElement(
                  "td",
                  null,
                  p.Median_per_cell_frit
                ),
                React.createElement(
                  "td",
                  null,
                  p.Median_duplication_rate
                )
              );
            })
          )
        )
      );
    }
  }]);

  return Table;
}(React.Component);

/*
** Main page layout.
*/


function ExperimentPage(props) {
  return React.createElement(
    "span",
    null,
    React.createElement(Header, { run_name: props.run_name }),
    React.createElement(
      "div",
      { className: "container-fluid" },
      React.createElement(
        "div",
        { className: "row" },
        React.createElement(
          "nav",
          { className: "col-md-2 d-none d-md-block bg-light sidebar" },
          React.createElement(
            "div",
            { className: "sidebar-sticky" },
            React.createElement(
              "div",
              { className: "nav flex-column nav-pills", id: "v-pills-tab", role: "tablist", "aria-orientation": "vertical" },
              React.createElement(
                "a",
                { className: "nav-link active", id: "summary-tab", "data-toggle": "pill", href: "#summary", role: "tab", "aria-controls": "summary", "aria-selected": "true" },
                "Summary Table"
              ),
              props.samp_pills
            )
          )
        ),
        React.createElement(
          "main",
          { role: "main", className: "col-md-9 ml-sm-auto col-lg-10 px-4", style: { paddingTop: "15px" } },
          React.createElement(
            "div",
            { className: "tab-content", id: "nav-tabContent" },
            React.createElement(Table, { data: tableData }),
            props.samp_list
          )
        )
      )
    )
  );
}

var sampList = run_data.sample_list.map(function (samp) {
  return React.createElement(Sample, { key: samp, sample_id: samp });
});

var sampPills = run_data.sample_list.map(function (samp) {
  return React.createElement(SamplePill, { key: samp, sample_id: samp });
});

ReactDOM.render(React.createElement(ExperimentPage, { samp_list: sampList, samp_pills: sampPills, run_name: run_data.run_name }), document.getElementById('exp_page'));