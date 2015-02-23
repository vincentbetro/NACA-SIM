// 
// JAVASCRIPT Functions for the JSE
//

// HTTPRequest variables
var xhs = false;
var xhs2 = false;
var xhs3 = false;

// Global variables for toggleAbstract()
var abstractStatus = new Array();
var abstractIDnum;

// Global variables for toggleRequestBox()
var requestStatus = new Array();
var requestIDnum;

// Global variables for processRequest()
var requestSubmitLock = "no";

// Global variables for setLinks()
var totalResults = 0;
var prev = -10;
var nextOne = 10;

// Preload image(s)
var loadbg = new Image();
loadbg.src = "bigloadicon.gif";


// Cookie read and write functions copied from http://www.quirksmode.org/js/cookies.html

function createCookie(name,value) {
	var date = new Date();
	// Set to expire in 12 hours
	date.setTime(date.getTime()+(12*60*60*1000));
	var expires = "; expires="+date.toGMTString();
	
	document.cookie = name+"="+value+expires+"; path=/";
}

function readCookie(name) {
	var nameEQ = name + "=";
	var ca = document.cookie.split(';');
	for(var i=0;i < ca.length;i++) {
		var c = ca[i];
		while (c.charAt(0)==' ') c = c.substring(1,c.length);
		if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length,c.length);
	}
	return null;
}

// restoreResultsFromCookie restores the contents of the searchResults DIV to the results previously saved

function restoreResultsFromCookie()
{
	if (readCookie("searchString"))
	{
		getSearchResults(1, 'all'); // Force getSearchResults to revert back to stored cookie data 
	}
}

// toggleResults toggles the visibility of the search results and search information

function toggleResults()
{
	var currBlockID = "results";
	if (document.getElementById(currBlockID).style.display == "block")
	{
		document.getElementById(currBlockID).style.display = "none";
		document.getElementById("searchNav").style.display = "none";
		document.getElementById("toggleResults").innerHTML = "<a href='javascript:toggleResults()'>Show</a>";
	}
	else
	{
		document.getElementById(currBlockID).style.display = "block";
		document.getElementById("searchNav").style.display = "block";
		document.getElementById("toggleResults").innerHTML = "<a href='javascript:toggleResults()'>Hide</a>";
	}
}


// entersSearch is called when "Enter" is pressed while the search box is focused.
// It responds by calling getSearchResults

function enterSearch(e)
{
	var keynum;

	if(window.event) // IE
  	{
  		keynum = e.keyCode;
  	}
	else if(e.which) // Netscape/Firefox/Opera
	{
  		keynum = e.which;
	}
	
	if (keynum == 13)
	{
		getSearchResults();
	}
}


// resetLimits resets the search result display limits

function resetLimits()
{
	prev = -10;
	nextOne = 10;
}


// getSearchResults provides AJAX functionality to the search function by sending the query and other information
// to search.php and calling getSearchOutput once the response from search.php is received. getSearchResults also
// stores previous search data in a cookie for retreival on reload

function getSearchResults(newOffset, useCookie)
{
	var encodedSearch = document.getElementById('searchString').value;
	encodedSearch = escape(encodedSearch);
	
	var journal = document.getElementById('journal').value;
	var queryType = "";
	
	// if useCookie is not set, use submitted search string and store new cookies
	if (!useCookie)
	{	
		if (document.getElementById("scope1").checked)
		{
			queryType = "both";
		}
		else if (document.getElementById("scope2").checked)
		{
			queryType = "title";
		}
		else
		{
			queryType = "fulltext";
		}
		
		createCookie("searchString", encodedSearch);
		createCookie("journal", journal);
		createCookie("queryType", queryType);
		createCookie("newOffset", newOffset);
	}
	// if useCookie is set to "searchString," just use the stored search string
	if (useCookie == "searchString")
	{
		encodedSearch = readCookie("searchString");
		journal = readCookie("journal");
		queryType = readCookie("queryType");
		createCookie("newOffset", newOffset);
	}
	// if useCookie is set to "all", use all data found in cookies
	if (useCookie == "all")
	{
		encodedSearch = readCookie("searchString");
		journal = readCookie("journal");
		queryType = readCookie("queryType");
		newOffset = readCookie("newOffset");
	}
	
	if (!(document.getElementById('toggleResults')))
	{
		document.getElementById("loadingImg").innerHTML = "<img src='loading.gif' alt='loading...' \>";
	}
	else
	{
		document.body.background = "bigloadicon.gif";
		document.getElementById("searchResults").style.filter="alpha(opacity=30)";
		document.getElementById("searchResults").style.MozOpacity=0.3;
	}	
	
	if (window.XMLHttpRequest)
	{
		xhs2 = new XMLHttpRequest();
	}
	else if (window.ActiveXObject)
	{
		try
		{
			xhs2 = new ActiveXObject("Microsoft.XMLHTTP");
		}
		catch (e)
		{}
	}

	if (xhs2)
	{
		
		// Set new limits
		nextOne = parseInt(newOffset) + 10;
		prev = parseInt(newOffset) - 10;

		var url = "search.php?searchString=" + encodedSearch + "&offset=" + newOffset + "&qType=" + queryType + "&journal=" + journal;
		xhs2.onreadystatechange = getSearchOutput;
		xhs2.open("GET", url, true);
		xhs2.send(null);
	}
	else
	{
		alert("Error 1!");
	}
}

// getSearchOutput sends the response from the HTTP search request to the "searchResults" HTML entity

function getSearchOutput()
{	
	if (xhs2.readyState == 4)
	{
		if (xhs2.status == 200)
		{
			document.getElementById('searchResults').innerHTML = xhs2.responseText;
			document.getElementById('loadingImg').innerHTML = "";
			document.body.background = "";
			document.getElementById("searchResults").style.filter="alpha(opacity=100)";
			document.getElementById("searchResults").style.MozOpacity="1.0";
			setLinks();
		}
		else
		{
			alert("Error 2!");
		}
	}
	else
	{
		// Do nothing
	}
}


// setLinks creates and displays the numerical and directional links for navigating the search results

function setLinks()
{
	var nextBlockID = "results";
	
	totalResults = parseInt(document.getElementById("resultNum").innerHTML);
	
	
	
	var totalUpToTen = 0; 
	if (totalResults % 10 == 0)
	{
		totalUpToTen = totalResults;
	}
	else
	{
		totalUpToTen = totalResults + 10 - (totalResults % 10);
	}		

	// Create the page bar

	var pageBar = "";
	var availableSlots = 21;

	if (totalResults > 10)
	{
		var curr = nextOne - 10;
		var begin = (prev - 80)/10;
		var end = (nextOne + 100)/10;
	
		var max = 21;
		if (totalUpToTen/10 < max)
		{
			max = totalUpToTen/10;
		}
	
		if (begin < 1)
		{
			begin = 1;
			end = max;
		}

		
		if (((totalUpToTen)/10) < end)
		{
			begin = (totalUpToTen)/10 - max + 1;
			end = (totalUpToTen)/10;
		}
		
		for (var i=begin; i <= end; i++)
		{
			if ((i-1)*10 < nextOne && (i-1)*10 > prev)
			{
				pageBar += " <span style='color: darkred'>" + i + "</span>";
			}
			else
			{
				pageBar += " <a href=\"javascript:getSearchResults(" + (i-1)*10 + ", 'searchString')\">" + i + "</a>";
			}
		}
	}
	
	var rangeEnd = 0;
	var rangeBegin = 0;
	if (nextOne < totalUpToTen && document.getElementById(nextBlockID))
	{
		rangeEnd = nextOne;
	}
	else
	{
		rangeEnd = totalResults;
	}
	if (prev < 0)
	{
		rangeBegin = 1;
	}
	else
	{
		rangeBegin = curr+1;
	}
	
	if (totalResults <= 10)
	{
		document.getElementById('searchNav').innerHTML = "";
		document.getElementById('searchNavBottom').innerHTML = "";
	}
	else if (prev < 0 && document.getElementById(nextBlockID))
	{
		document.getElementById('searchNav').innerHTML = "<strong>Back : " + pageBar + " : <a href='javascript:getNextBlock()'>Next</a> &nbsp;&nbsp;(" + rangeBegin + " to " + rangeEnd + ")</strong>";
		document.getElementById('searchNavBottom').innerHTML = "<strong>Back : " + pageBar + " : <a href='javascript:getNextBlock()'>Next</a> &nbsp;&nbsp;(" + rangeBegin + " to " + rangeEnd + ")</strong>";
	}
	else if (prev >= 0 && nextOne < totalUpToTen && document.getElementById(nextBlockID))
	{
		document.getElementById('searchNav').innerHTML = "<strong><a href='javascript:getPrevBlock()'>Back</a> : " + pageBar + " : <a href='javascript:getNextBlock()'>Next</a> &nbsp;&nbsp;(" + rangeBegin + " to " + rangeEnd + ")</strong>";
		document.getElementById('searchNavBottom').innerHTML = "<strong><a href='javascript:getPrevBlock()'>Back</a> : " + pageBar + " : <a href='javascript:getNextBlock()'>Next</a> &nbsp;&nbsp;(" + rangeBegin + " to " + rangeEnd + ")</strong>";
	}
	else if (prev >= 0 && nextOne >= totalUpToTen && document.getElementById(nextBlockID))
	{
		document.getElementById('searchNav').innerHTML = "<strong><a href='javascript:getPrevBlock()'>Back</a> : " + pageBar + " : Next &nbsp;&nbsp;(" + rangeBegin + " to " + rangeEnd + ")</strong>";
		document.getElementById('searchNavBottom').innerHTML = "<strong><a href='javascript:getPrevBlock()'>Back</a> : " + pageBar + " : Next &nbsp;&nbsp;(" + rangeBegin + " to " + rangeEnd + ")</strong>";
	}
	else
	{
		if (document.getElementById('searchNav'))
		{
			document.getElementById('searchNav').innerHTML = "";
			document.getElementById('searchNavBottom').innerHTML = "";
		}
	}

}


// getPrevBlock is called by clicking "Back" and simply makes a call to getSearchResults requesting 
// the previous 10 results

function getPrevBlock()
{
	getSearchResults(prev, 'searchString');
}


// getNextBlock is called by clicking "Next" and simply makes a call to getSearchResults requesting 
// the next 10 results

function getNextBlock()
{
	
	if (document.getElementById("results"))
	{
		getSearchResults(nextOne, 'searchString');
	}
}


// hideList is useless and outdated. It should probably be deleted

function hideList(listID)
{
	if (document.getElementById(listID).style.visibility != "hidden")
	{
		document.getElementById(listID).style.visibility = "hidden";
		document.getElementById(listID).style.display = "none";
	}
	else
	{
		document.getElementById(listID).style.visibility = "visible";
		document.getElementById(listID).style.display = "inline";
	}
}


// expandLoginInfo simply makes the login submission form visible

function expandLoginInfo()
{
	document.getElementById('addJournalInfo').innerHTML = "<br />\n <table width='550' border='0'>\n <tr> <td><b>Journal Provider URL:</b></td>\n <td><input type='text' size='50' id='provider' name='provider' value='http://'></td>\n </tr><tr>\n <td><b>Username:</b></td>\n <td><input type='text' id='username' name='username'></td>\n </tr><tr>\n <td><b>Password:</b></td>\n <td><input type='password' id='password' name='password'></td>\n </tr>\n </table>\n <br />\n Your password will be saved in a permission-limited encrypted file.<br /><br /> <input type='submit' value='Submit Info'><br /><br /> <a href='javascript:hideLoginInfo()'><b>&lt;&lt; Hide</b></a>";
}


// hideLoginInfo simply makes the login submission form hidden

function hideLoginInfo()
{
	document.getElementById('addJournalInfo').innerHTML = "<b><a href='javascript:expandLoginInfo()'>Submit Login Info &gt;&gt;</a></b>"
}


// validateLoginInfo validates the login submission form and returns a boolean value determining whether
// or not the form can be submitted

function validateLoginInfo()
{
	var provider = document.getElementById('provider').value;
	var username = document.getElementById('username').value;
	var password = document.getElementById('password').value;
	
	document.getElementById('provider').style.background = "#ffffff";
	document.getElementById('username').style.background = "#ffffff";
	document.getElementById('password').style.background = "#ffffff";
	
	var valid = "yes";
	
	if (provider == "http://")
	{
		document.getElementById('provider').style.background = "red";
		valid = "no";
	}
	
	if (username == "")
	{
		document.getElementById('username').style.background = "red";
		valid = "no";
	}
	
	if (password == "")
	{
		document.getElementById('password').style.background = "red";
		valid = "no";
	}
	
	if (valid == "yes")
	{
		return true;
	}
	else
	{
		return false;
	}
}


// toggleAbstract expands and closes the abstract for a particular article listing

function toggleAbstract(aid)
{
	var abstractID = aid + '-abstract';
	var aImageID = aid + '-aimage';

	// Setting status to state 2 means that the abstract has been fetched and displayed
	// Setting status to state 1 means that the abstract has been fetched but is hidden
	// Default status of implied state 0 means no data has yet been fetched or displayed

	if (abstractStatus[abstractID] == 2) 
	{
		new Effect.BlindUp(abstractID);
		abstractStatus[abstractID] = 1;
		document.getElementById(aImageID).src = "icons/aicon.png";
	}
	else if (abstractStatus[abstractID] == 1)
	{
		new Effect.BlindDown(abstractID);
		abstractStatus[abstractID] = 2;
		document.getElementById(aImageID).src = "icons/aicon-selected.png";
	}
	else
	{
		new Effect.BlindDown(abstractID);
		getAbstract(aid);
		abstractStatus[abstractID] = 2;
		document.getElementById(aImageID).src = "icons/aicon-selected.png";
	}
}

// getAbstract provides AJAX functionality to the "Abstract" button, allowing
// the database to be updated without refreshing the page 

function getAbstract(aid)
{
	
	if (window.XMLHttpRequest)
	{
		xhs3 = new XMLHttpRequest();
	}
	else if (window.ActiveXObject)
	{
		try
		{
			xhs3 = new ActiveXObject("Microsoft.XMLHTTP");
		}
		catch (e)
		{}
	}

	if (xhs3)
	{
		abstractIDnum = aid;

		var url = "getAbstract.pl?id=" + aid;
		xhs3.onreadystatechange = returnAbstract;
		xhs3.open("GET", url, true);
		xhs3.send(null);
	}
	else
	{
		alert("Error 1!");
	}
}

// confirmRequest changes the appearance of the "request loan" button after the 
// request has been made

function returnAbstract()
{	
	if (xhs3.readyState == 4)
	{
		if (xhs3.status == 200)
		{
			var spanID = abstractIDnum + "-abstract";
			document.getElementById(spanID).innerHTML = xhs3.responseText;
		}
		else
		{
			alert("Error 2!");
		}
	}
	else
	{
		// Do nothing
	}
}


// toggleRequest expands and closes the request box for a particular article listing

function toggleRequestBox(rid)
{
	var requestID = rid + '-request';
	var rImageID = rid + '-rimage';

	// Setting status to state 2 means that the request box has been fetched and displayed
	// Setting status to state 1 means that the request box has been fetched but is hidden
	// Default status of implied state 0 means no data has yet been fetched or displayed

	if (requestStatus[requestID] == 2) 
	{
		new Effect.BlindUp(requestID);
		requestStatus[requestID] = 1;
		document.getElementById(rImageID).src = "icons/requesticon.png";
	}
	else if (requestStatus[requestID] == 1)
	{
		new Effect.BlindDown(requestID);
		requestStatus[requestID] = 2;
		document.getElementById(rImageID).src = "icons/requesticon-selected.png";
	}
	else
	{
		getRequestBox(rid);
		new Effect.BlindDown(requestID);
		requestStatus[requestID] = 2;
		document.getElementById(rImageID).src = "icons/requesticon-selected.png";
	}
}

// getRequestBox generates the contents of the request box (without AJAX)

function getRequestBox(rid)
{
	var boxContent = "<br/><strong><span id=\"" + rid + "-rstatus\">Notice: Physical copies may already be present in the Lupton Library</span></strong>";
	boxContent += "<br/><em>Inter-library loans require a UTC account <a href=\"http://illiad.lib.utc.edu/illiad.dll\">set up with ILLiad</a></em>:<br/>";
	boxContent += "<strong>UTCID:</strong> <input type=\"text\" id=\"utcid\" size=\"6\" onkeypress=\"enterRequest(event," + rid + ")\"> <strong>Password:</strong> <input type=\"password\" id=\"utcpass\" size=\"12\" onkeypress=\"enterRequest(event," + rid + ")\"> ";
	boxContent += "<input type=\"button\" value=\"Send Request\" onclick=\"processRequest(" + rid + ")\"><br/><br/>";

	var spanID = rid + "-request";
	document.getElementById(spanID).innerHTML = boxContent;
}

// processRequest sends the request login info and data to a submission script which checks for valid login
// before passing the data to ILLiad (uses AJAX)

function processRequest(rid)
{
	if (requestSubmitLock == "no")
	{	
		if (window.XMLHttpRequest)
		{
			xhs3 = new XMLHttpRequest();
		}
		else if (window.ActiveXObject)
		{
			try
			{
				xhs3 = new ActiveXObject("Microsoft.XMLHTTP");
			}
			catch (e)
			{}
		}
	
		if (xhs3)
		{
			requestSubmitLock = "yes";	

			requestIDnum = rid;
			var user = escape(document.getElementById("utcid").value);
			var pass = escape(document.getElementById("utcpass").value);
	
			var url = "processRequest.pl?id=" + rid + "&user=" + user + "&pass=" + pass;
			xhs3.onreadystatechange = returnRequestStatus;
			xhs3.open("GET", url, true);
			xhs3.send(null);
		}
		else
		{
			alert("Error 1!");
		}
	}
}

// doNothing does nothing
function doNothing()
{}

// returnRequestStatus returns the status of a processed request

function returnRequestStatus()
{	
	if (xhs3.readyState == 4)
	{
		if (xhs3.status == 200)
		{
			var requestStatus = xhs3.responseText;
			var statusSpan = requestIDnum + '-rstatus';

			if (requestStatus == "Wrong Password!")
			{
				requestSubmitLock = "no";
				document.getElementById(statusSpan).innerHTML = "<span style='color:red'>Wrong Password!</span>";
			}
			else
			{
				document.getElementById(statusSpan).innerHTML = "Request Sent!";

				var toggleString = "toggleRequestBox(" + requestIDnum + ")";
				setTimeout(toggleString, 1000);
				//toggleRequestBox(requestIDnum);
				requestSubmitLock = "no";
			}
		}
		else
		{
			alert("Error 2!");
		}
	}
	else
	{
		// Do nothing
	}
}

// entersSearch is called when "Enter" is pressed while a request text field is focused.
// It responds by calling processRequest()

function enterRequest(e, rid)
{
	var keynum;

	if(window.event) // IE
  	{
  		keynum = e.keyCode;
  	}
	else if(e.which) // Netscape/Firefox/Opera
	{
  		keynum = e.which;
	}
	
	if (keynum == 13)
	{
		processRequest(rid);
	}
}


