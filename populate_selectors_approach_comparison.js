// speed - populate selectors
// Outcomes:
// 1. '<option value=...' is a bit faster than new Option
// 2. innerHTML is x1.5 better than JQuery append
// 3. for of vs native for (length chached, var used) vs forEach - no major difference
// 4. pJqueryMap is second slowest after appending on each iteration

let arr = [...Array(99999).keys()];
arr = arr.map(i => i.toString());
let ele = instance.domainsSelectEle;

let functionArr = [];

functionArr[0] = function pJqueryMap(ele, optionsArr) {
				ele.options.length = 0;
				var start = new Date();
				
				$(ele).append($.map(optionsArr, (item) => `<option value="${item}">${item}</option>`));
				
				var time = new Date() - start;
				console.log(time);
				
				$(ele).selectpicker('refresh');
				return time;
			}
			
functionArr[1] = function pForOfAppend(ele, optionsArr) {
				ele.options.length = 0;
				var start = new Date();
				
				for (const item of optionsArr) {
					$(ele).append(new Option(item, item))
				}
				
				var time = new Date() - start;
				console.log(time);
				
				$(ele).selectpicker('refresh');
				return time;
			}
			
functionArr[2] = function pForOf(ele, optionsArr) {
				ele.options.length = 0;
				var start = new Date();
				
				let html = '';
				for (const item of optionsArr) {
					html = html + `<option value="${item}">${item}</option>`;
				}
				$(ele).append(html);
				
				var time = new Date() - start;
				console.log(time);
				
				$(ele).selectpicker('refresh');
				return time;
			}

functionArr[3] = function pForTraditionalLChachedInnerHTML(ele, optionsArr) {
				ele.options.length = 0;
				var start = new Date();
				
				let html = '';
				// superoptimized
				for (var i = 0, len = optionsArr.length; i < len; ++i) {
					html = html + `<option value="${optionsArr[i]}">${optionsArr[i]}</option>`;
				}
				ele.innerHTML = html;
				
				var time = new Date() - start;
				console.log(time);
				
				$(ele).selectpicker('refresh');
				return time;
			}		
		
functionArr[4] = function pForTraditionalLChached(ele, optionsArr) {
				ele.options.length = 0;
				var start = new Date();
				
				let html = '';
				// superoptimized
				for (var i = 0, len = optionsArr.length; i < len; ++i) {
					html = html + `<option value="${optionsArr[i]}">${optionsArr[i]}</option>`;
				}
				$(ele).append(html);
				
				var time = new Date() - start;
				console.log(time);
				
				$(ele).selectpicker('refresh');
				return time;
			}

functionArr[5] = function pForTraditionalLChachedOptionConstructor(ele, optionsArr) {
				ele.options.length = 0;
				var start = new Date();
				
				let html = [];
				// superoptimized
				for (var i = 0, len = optionsArr.length; i < len; ++i) {
					html.push(new Option(optionsArr[i], optionsArr[i]));
				}
				$(ele).append(html);
				
				var time = new Date() - start;
				console.log(time);
				
				$(ele).selectpicker('refresh');
				return time;
			}

functionArr[6] = function pForEach(ele, optionsArr) {
				ele.options.length = 0;
				var start = new Date();
				
				let html = '';
				optionsArr.forEach((item) => {
					html = html + `<option value="${item}">${item}</option>`;
				})
				
				$(ele).append(html);
				
				var time = new Date() - start;
				console.log(time);
				
				$(ele).selectpicker('refresh');
				return time;
			}

let timeArr = [];	