/* File : tables.js
 *
 * Contains onchange() functions for the pairwise comparison table and risk of
 * bias tables. The function are bound by formatting functions in `util.R`.
 *
 * See: `FormatColumnStudyEffects()`, `FormatPairwiseTableDropdown()` and
 * `FormatColumnContrEval()`
 */
function selectPairwiseTableItem(id) {
  var x = document.getElementById(id).value;
  console.log("you selected",id);
  console.log("value",x);
  var out = {id:id, value:x};
  Shiny.setInputValue("pairwiseSelect", out);
}

function changeContrEval(id) {
  var x = document.getElementById(id).value;
  console.log("you selected",id);
  console.log("value",x);
  var out = {id:id, value:x};
  Shiny.setInputValue("contrEval", out);
}

function changeStudyEffects(id) {
  var x = document.getElementById(id).value;
  var out = {id:id, value:x};
  Shiny.setInputValue("studyEffects", out);
  console.log("Set value of " + id + " to " + x);
}

function changeOverallRob(id) {
  var x = document.getElementById(id).value;
  var out = { id: id, value: x };
  Shiny.setInputValue("overallRob", out);
}
